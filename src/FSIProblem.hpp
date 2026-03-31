#ifndef FSI_PROBLEM_HPP
#define FSI_PROBLEM_HPP

#include <deal.II/base/conditional_ostream.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/distributed/fully_distributed_tria.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_renumbering.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/fe/fe_nothing.h>
#include <deal.II/fe/fe_simplex_p.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/fe_values_extractors.h>
#include <deal.II/fe/mapping_fe.h>
#include <deal.II/grid/grid_in.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/hp/fe_collection.h>
#include <deal.II/hp/fe_values.h>
#include <deal.II/hp/mapping_collection.h>
#include <deal.II/hp/q_collection.h>
#include <deal.II/lac/affine_constraints.h>
#include <deal.II/lac/block_sparse_matrix.h>
#include <deal.II/lac/block_vector.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/solver_gmres.h>
#include <deal.II/lac/sparse_direct.h>
#include <deal.II/lac/trilinos_block_sparse_matrix.h>
#include <deal.II/lac/trilinos_parallel_block_vector.h>
#include <deal.II/lac/trilinos_precondition.h>
#include <deal.II/lac/trilinos_solver.h>
#include <deal.II/lac/trilinos_sparse_matrix.h>
#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/numerics/vector_tools.h>

#include <fstream>
#include <iostream>

using namespace dealii;

/**
 * Class managing the differential problem.
 */
class FSIProblem {
   public:
    // Physical dimension (1D, 2D, 3D)
    static constexpr unsigned int dim = 2;

    // Enum for boundary Ids
    enum BoundaryIds : types::boundary_id {
        Inflow = 14,
        Outflow = 15,
        BottomWall = 16,
        LeftWallSolid = 17,
        LeftWallFluid = 18,
        RightWallSolid = 19,
        RightWallFluid = 20
    };

    // Enum for material Ids
    enum MaterialIds : types::material_id { Fluid = 12, Solid = 13 };

    class InletVelocity : public Function<dim> {
       public:
        InletVelocity() : Function<dim>(dim + 1) {}

        virtual void vector_value(const Point<dim>& p,
                                  Vector<double>& values) const override {
            values[0] = 0.0;                                       // u_x
            values[1] = (std::sin(M_PI * (2.0 * p[0]))) * (-3.0);  // u_y

            values[2] = 0.0;
        }

        virtual double value(const Point<dim>& p,
                             const unsigned int component = 0) const override {
            if (component == 1)
                return std::sin(M_PI * (2.0 * p[0])) * (-3.0);
            else
                return 0.0;
        }
    };

    // Constructor.
    FSIProblem(const std::string& mesh_file_name_,
               const unsigned int& degree_velocity_,
               const unsigned int& degree_pressure_,
               const unsigned int& degree_displacement_)
        : mesh_file_name(mesh_file_name_),
          degree_velocity(degree_velocity_),
          degree_pressure(degree_pressure_),
          degree_displacement(degree_displacement_),
          // for parallelization
          mpi_size(Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD)),
          mpi_rank(Utilities::MPI::this_mpi_process(MPI_COMM_WORLD)),
          pcout(std::cout, mpi_rank == 0),  // only rank=0 can write
          mesh(MPI_COMM_WORLD) {}
    // Initialization.
    void setup();

    // System assembly.
    void assemble();

    // System solution.
    void solve();

    // Output.
    void output();

   protected:
    // Problem definition.
    // ///////////////////////////////////////////////////////

    // Kinematic viscosity [m2/s].
    const double nu = 1.0;

    // Inlet velocity.
    InletVelocity inlet_velocity;

    // Lamé parameters.
    const double mu = 1.0;
    const double lambda = 1.0;

    // Discretization.
    // ///////////////////////////////////////////////////////////

    // Name of the mesh.
    const std::string mesh_file_name;

    // Polynomial degree used for velocity.
    const unsigned int degree_velocity;

    // Polynomial degree used for pressure.
    const unsigned int degree_pressure;

    // Polynomial degree used for displacement.
    const unsigned int degree_displacement;

    // the order of declaration of these  4variables is important
    const unsigned int mpi_size;
    const unsigned int mpi_rank;

    ConditionalOStream pcout;

    // Mesh.
    parallel::fullydistributed::Triangulation<dim> mesh;

    // hp collections
    hp::FECollection<dim> fe_collection;
    hp::QCollection<dim> quadrature_collection;
    hp::QCollection<dim - 1> face_quadrature_collection;
    hp::MappingCollection<dim> mapping_collection;

    // Constraints
    AffineConstraints<double> constraints;

    // fe
    std::unique_ptr<FESystem<dim>> fe_stokes;
    std::unique_ptr<FESystem<dim>> fe_elasticity;

    // Sparsity
    SparsityPattern sparsity;

    // DoF handler.
    DoFHandler<dim> dof_handler;

    // DoFs owned by current process.
    IndexSet locally_owned_dofs;

    // DoFs relevant to the current process (including ghost DoFs).
    IndexSet locally_relevant_dofs;

    // System matrix.
    TrilinosWrappers::SparseMatrix system_matrix;

    // Right-hand side vector in the linear system.
    TrilinosWrappers::MPI::Vector system_rhs;

    // System solution (without ghost elements).
    TrilinosWrappers::MPI::Vector solution_owned;

    // System solution (including ghost elements).
    TrilinosWrappers::MPI::Vector solution;
};

#endif
