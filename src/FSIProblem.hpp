#ifndef FSI_PROBLEM_HPP
#define FSI_PROBLEM_HPP

#include <deal.II/base/conditional_ostream.h>
#include <deal.II/base/quadrature_lib.h>

#include <deal.II/distributed/fully_distributed_tria.h>

#include <deal.II/dofs/dof_renumbering.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/dofs/dof_handler.h> 
#include <deal.II/lac/affine_constraints.h>

#include <deal.II/fe/fe_simplex_p.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/mapping_fe.h>
#include <deal.II/fe/fe_values_extractors.h>
#include <deal.II/fe/fe_nothing.h>
#include <deal.II/hp/fe_collection.h>
#include <deal.II/hp/q_collection.h>
#include <deal.II/hp/mapping_collection.h>
#include <deal.II/hp/fe_values.h>

#include <deal.II/grid/grid_in.h>
#include <deal.II/grid/grid_tools.h>

#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/solver_gmres.h>
#include <deal.II/lac/sparse_direct.h>
#include <deal.II/lac/trilinos_block_sparse_matrix.h>
#include <deal.II/lac/trilinos_parallel_block_vector.h>
#include <deal.II/lac/trilinos_precondition.h>
#include <deal.II/lac/trilinos_sparse_matrix.h>
#include <deal.II/lac/block_sparse_matrix.h>
#include <deal.II/lac/block_vector.h>

#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/numerics/vector_tools.h>

#include <fstream>
#include <iostream>

using namespace dealii;

/**
 * Class managing the differential problem.
 */
class FSIProblem
{
public:
  // Physical dimension (1D, 2D, 3D)
  static constexpr unsigned int dim = 2;

  // Enum for boundary Ids
  enum BoundaryIds : types::boundary_id
  {
    Inflow = 14,
    Outflow = 15,
    BottomWall = 16,
    LeftWallSolid = 17,
    LeftWallFluid = 18,
    RightWallSolid = 19,
    RightWallFluid = 20
  };

  // Enum for material Ids
  enum MaterialIds : types::material_id
  {
    Fluid = 12,
    Solid = 13
  };

  class InletVelocity : public Function<dim>
  {
  public:
    InletVelocity()
      : Function<dim>(dim + 1)
    {}

    virtual void
    vector_value(const Point<dim> &p, Vector<double> &values) const override
    {
      values[0] = 0.0; // u_x
      values[1] = (std::sin(M_PI * (2.0 *p[0])))*(-3.0); // u_y

      values[2] = 0.0;
    }

    virtual double
    value(const Point<dim> &p, const unsigned int component = 0) const override
    {
      if (component == 1)
        return std::sin(M_PI * (2.0 * p[0]))*(-3.0);
      else
        return 0.0;
    }
  };


  class PreconditionIdentity
  {
  public:
    // Application of the preconditioner: we just copy the input vector (src)
    // into the output vector (dst).
    void
    vmult(Vector<double>       &dst,
          const Vector<double> &src) const
    {
      dst = src;
    }

  protected:
  };

  /* Block-diagonal preconditioner.
  class PreconditionBlockDiagonal
  {
  public:
    // Initialize the preconditioner, given the velocity stiffness matrix, the
    // pressure mass matrix.
    void
    initialize(const TrilinosWrappers::SparseMatrix &velocity_stiffness_,
               const TrilinosWrappers::SparseMatrix &pressure_mass_,
               const TrilinosWrappers::SparseMatrix &displacement_stiffness_)
    {
      velocity_stiffness = &velocity_stiffness_;
      pressure_mass      = &pressure_mass_;
      displacement_stiffness = &displacement_stiffness_;

      preconditioner_velocity.initialize(velocity_stiffness_);
      preconditioner_pressure.initialize(pressure_mass_);
      preconditioner_displacement.initialize(displacement_stiffness_);
    }

    // Application of the preconditioner.
    void
    vmult(TrilinosWrappers::MPI::BlockVector       &dst,
          const TrilinosWrappers::MPI::BlockVector &src) const
    {
      SolverControl solver_control_velocity(1000, 1e-2 * src.block(0).l2_norm());
      SolverCG<TrilinosWrappers::MPI::Vector> solver_cg_velocity(
        solver_control_velocity);
      solver_cg_velocity.solve(*velocity_stiffness,
                               dst.block(0),
                               src.block(0),
                               preconditioner_velocity);

      SolverControl                           solver_control_pressure(1000,
                                            1e-2 * src.block(1).l2_norm());
      SolverCG<TrilinosWrappers::MPI::Vector> solver_cg_pressure(
        solver_control_pressure);
      solver_cg_pressure.solve(*pressure_mass,
                               dst.block(1),
                               src.block(1),
                               preconditioner_pressure);

      SolverControl                           solver_control_displacement(1000,
                                            1e-2 * src.block(2).l2_norm());
      SolverCG<TrilinosWrappers::MPI::Vector> solver_cg_displacement(
        solver_control_displacement);
      solver_cg_displacement.solve(*displacement_stiffness,
                               dst.block(2),
                               src.block(2),
                               preconditioner_displacement);
    }

  protected:
    // Velocity stiffness matrix.
    const TrilinosWrappers::SparseMatrix *velocity_stiffness;

    // Preconditioner used for the velocity block.
    TrilinosWrappers::PreconditionILU preconditioner_velocity;

    // Pressure mass matrix.
    const TrilinosWrappers::SparseMatrix *pressure_mass;

    // Preconditioner used for the pressure block.
    TrilinosWrappers::PreconditionILU preconditioner_pressure;

    // displacement stiffness matrix.
    const TrilinosWrappers::SparseMatrix *displacement_stiffness;

    // Preconditioner used for the displacement block.
    TrilinosWrappers::PreconditionILU preconditioner_displacement;

  };
*/
  // Constructor.
  FSIProblem(const std::string  &mesh_file_name_,
                    const unsigned int &degree_velocity_,
                    const unsigned int &degree_pressure_,
                    const unsigned int &degree_displacement_)
    : mpi_size(Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD))
    , mpi_rank(Utilities::MPI::this_mpi_process(MPI_COMM_WORLD))
    , pcout(std::cout, mpi_rank == 0)
    , mesh_file_name(mesh_file_name_)
    , degree_velocity(degree_velocity_)
    , degree_pressure(degree_pressure_)
    , degree_displacement(degree_displacement_)
    , mesh()
  {}

  // Initialization.
  void
  setup();

  // System assembly.
  void
  assemble();

  // System solution.
  void
  solve();

  // Output.
  void
  output();

protected:

  // MPI parallel. /////////////////////////////////////////////////////////////

  // Number of MPI processes.
  const unsigned int mpi_size;

  // This MPI process.
  const unsigned int mpi_rank;

  // Parallel output stream.
  ConditionalOStream pcout;

  // Problem definition. ///////////////////////////////////////////////////////

  // Kinematic viscosity [m2/s].
  const double nu = 1.0;

  // Outlet pressure [Pa].
  //const double p_out = 0.0;

  // Inlet velocity.
  InletVelocity inlet_velocity;

  // Lamé parameters.
  const double mu     = 1.0;
  const double lambda = 1.0;


  // Discretization. ///////////////////////////////////////////////////////////

  // Name of the mesh.
  const std::string mesh_file_name;

  // Polynomial degree used for velocity.
  const unsigned int degree_velocity;

  // Polynomial degree used for pressure.
  const unsigned int degree_pressure;

  // Polynomial degree used for displacement.
  const unsigned int degree_displacement;

  // Mesh.
  Triangulation<dim> mesh;

  // hp collections
  hp::FECollection<dim> fe_collection;
  hp::QCollection<dim> quadrature_collection;
  hp::QCollection<dim-1> face_quadrature_collection;
  hp::MappingCollection<dim> mapping_collection;

  // Constraints
  AffineConstraints<double> constraints;

  // fe
  std::unique_ptr<FESystem<dim>> fe_stokes;
  std::unique_ptr<FESystem<dim>> fe_elasticity;

  // Sparsity
  SparsityPattern sparsity;
  SparsityPattern sparsity_pressure_mass;

  // DoF handler.
  DoFHandler<dim> dof_handler;

  // DoFs owned by current process.
  IndexSet locally_owned_dofs;

  // DoFs owned by current process in the velocity, pressure and displacement blocks.
  std::vector<IndexSet> block_owned_dofs;

  // DoFs relevant to the current process (including ghost DoFs).
  IndexSet locally_relevant_dofs;

  // DoFs relevant to current process in the velocity, pressure and displacement blocks.
  std::vector<IndexSet> block_relevant_dofs;

  // System matrix.
  SparseMatrix<double> system_matrix;

  //SparseMatrix<double> pressure_mass;

  // Right-hand side vector in the linear system.
  Vector<double> system_rhs;

  // System solution (without ghost elements).
  Vector<double> solution_owned;

  // System solution (including ghost elements).
  Vector<double> solution;

};

#endif