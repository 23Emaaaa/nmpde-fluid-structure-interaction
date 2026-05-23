#ifndef FSI_PROBLEM_HPP
#define FSI_PROBLEM_HPP

#include <deal.II/base/conditional_ostream.h>
#include <deal.II/base/mpi.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/utilities.h>

#include <deal.II/distributed/shared_tria.h>

#include <deal.II/dofs/dof_renumbering.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/lac/affine_constraints.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/sparsity_tools.h>

#include <deal.II/fe/fe_simplex_p.h>
#include <deal.II/fe/fe_nothing.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/fe_values_extractors.h>
#include <deal.II/fe/mapping_fe.h>

#include <deal.II/hp/fe_collection.h>
#include <deal.II/hp/fe_values.h>
#include <deal.II/hp/mapping_collection.h>
#include <deal.II/hp/q_collection.h>

#include <deal.II/grid/grid_in.h>
#include <deal.II/grid/grid_tools.h>

#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/solver_gmres.h>
#include <deal.II/lac/trilinos_block_sparse_matrix.h>
#include <deal.II/lac/trilinos_parallel_block_vector.h>
#include <deal.II/lac/trilinos_precondition.h>
#include <deal.II/lac/trilinos_sparsity_pattern.h>
#include <deal.II/lac/trilinos_sparse_matrix.h>
#include <deal.II/lac/trilinos_vector.h>

#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/vector_tools.h>

#include <fstream>
#include <iostream>
#include <memory>
#include <string>
#include <vector>

using namespace dealii;

/**
 * Monolithic fluid-structure interaction problem on a simplicial mesh.
 */
class FSIProblem
{
public:
  // Physical dimension of the problem (1D, 2D or 3D).
  static constexpr unsigned int dim = 2;

  // Boundary identifiers expected from the mesh.
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

  // Material identifiers used to switch between fluid and solid cells.
  enum MaterialIds : types::material_id
  {
    Fluid = 12,
    Solid = 13
  };

  // Parabolic inlet profile prescribed on the fluid velocity block.
  class InletVelocity : public Function<dim>
  {
  public:
    InletVelocity()
        : Function<dim>(2 * dim + 1)
    {
    }

    virtual void vector_value(const Point<dim> &p, Vector<double> &values) const override
    {
      const double L = 0.5;
      const double Umax = 1.0;
      const double xi = p[0] / L;

      values = 0.0;
      values[0] = 0.0;                           // u_x
      values[1] = -4.0 * Umax * xi * (1.0 - xi); // u_y (downwards)
    }

    virtual double
    value(const Point<dim> &p, const unsigned int component = 0) const override
    {
      if (component != 1)
        return 0.0;

      const double L = 0.5;
      const double Umax = 1.0;
      const double xi = p[0] / L;

      return -4.0 * Umax * xi * (1.0 - xi);
    }
  };

  /**
   * Build the problem from the mesh path and polynomial degrees.
   */
  FSIProblem(const std::string &mesh_file_name_,
             const unsigned int &degree_velocity_,
             const unsigned int &degree_pressure_,
             const unsigned int &degree_displacement_)
      : mpi_communicator(MPI_COMM_WORLD),
        mpi_size(Utilities::MPI::n_mpi_processes(mpi_communicator)),
        mpi_rank(Utilities::MPI::this_mpi_process(mpi_communicator)),
        pcout(std::cout, mpi_rank == 0),
        mesh_file_name(mesh_file_name_),
        degree_velocity(degree_velocity_),
        degree_pressure(degree_pressure_),
        degree_displacement(degree_displacement_),
        mesh(mpi_communicator)
  {
  }

  // Read the mesh and initialize FE spaces, DoFs, constraints, and matrices.
  void setup();

  // Assemble the coupled fluid-solid linear system.
  void assemble();

  // Solve the linear system and store the distributed solution.
  void solve();

  // Write the solution fields in a visualization-friendly format.
  void output();

  // Write solver iterations and field statistics to a text report.
  void write_diagnostics_report(const std::string &output_file_name) const;

  // Gather the distributed solution on rank 0 and dump it to a text file.
  void write_solution_snapshot(const std::string &output_file_name) const;

protected:
  // MPI parallel.

  // MPI communicator.
  const MPI_Comm mpi_communicator;

  // Number of MPI processes.
  const unsigned int mpi_size;

  // This MPI process.
  const unsigned int mpi_rank;

  // Rank-0-only output stream.
  ConditionalOStream pcout;

  // Problem definition.

  // Fluid kinematic viscosity.
  const double nu = 1.0;

  // Prescribed inlet profile.
  InletVelocity inlet_velocity;

  // Solid Lamé parameters.
  const double mu = 1.0;
  const double lambda = 1.0;

  // Discretization.

  // Name of the mesh.
  const std::string mesh_file_name;

  // Polynomial degree used for velocity.
  const unsigned int degree_velocity;

  // Polynomial degree used for pressure.
  const unsigned int degree_pressure;

  // Polynomial degree used for displacement.
  const unsigned int degree_displacement;

  // Parallel shared mesh, which supports the simplicial Gmsh grid used.
  parallel::shared::Triangulation<dim> mesh;

  // hp collections used to switch between fluid and solid discretizations.
  hp::FECollection<dim> fe_collection;
  hp::QCollection<dim> quadrature_collection;
  hp::QCollection<dim - 1> face_quadrature_collection;
  hp::MappingCollection<dim> mapping_collection;

  // Global affine constraints.
  AffineConstraints<double> constraints;

  // FE systems for the fluid and solid subproblems.
  std::unique_ptr<FESystem<dim>> fe_stokes;
  std::unique_ptr<FESystem<dim>> fe_elasticity;

  // Block sparsity pattern of the coupled system.
  TrilinosWrappers::BlockSparsityPattern sparsity;

  // DoF handler.
  DoFHandler<dim> dof_handler;

  // DoFs owned by current process.
  IndexSet locally_owned_dofs;

  // DoFs relevant to the current process (including ghost DoFs).
  IndexSet locally_relevant_dofs;

  // DoFs owned by current process in the velocity, pressure and displacement blocks.
  std::vector<IndexSet> block_owned_dofs;

  // DoFs relevant to current process in the velocity, pressure and displacement blocks.
  std::vector<IndexSet> block_relevant_dofs;

  // Global DoF offsets for [u, p, d] after component-wise renumbering.
  std::vector<types::global_dof_index> block_offsets;

  // System matrix.
  TrilinosWrappers::BlockSparseMatrix system_matrix;

  // Pressure mass matrix used for the Stokes Schur complement preconditioner.
  TrilinosWrappers::SparseMatrix pressure_mass_matrix;

  // Right-hand side vector in the linear system.
  TrilinosWrappers::MPI::BlockVector system_rhs;

  // System solution (without ghost elements).
  TrilinosWrappers::MPI::BlockVector solution_owned;

  // System solution (including ghost elements).
  TrilinosWrappers::MPI::BlockVector solution;

  // Iteration counters from the last solve, used in diagnostics/reporting.
  unsigned int last_stokes_iterations = 0;
  unsigned int last_solid_iterations = 0;
};

#endif
