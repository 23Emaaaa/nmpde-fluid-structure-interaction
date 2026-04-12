#ifndef FSIPROBLEM_HPP
#define FSIPROBLEM_HPP

#include <deal.II/base/conditional_ostream.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/function.h>

#include <deal.II/distributed/fully_distributed_tria.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/dofs/dof_renumbering.h>

#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_simplex_p.h>
#include <deal.II/fe/fe_nothing.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/fe_values_extractors.h>
#include <deal.II/fe/mapping_fe.h>

#include <deal.II/grid/grid_in.h>
#include <deal.II/grid/grid_tools.h>

#include <deal.II/hp/fe_collection.h>
#include <deal.II/hp/fe_values.h>
#include <deal.II/hp/q_collection.h>
#include <deal.II/hp/mapping_collection.h>

#include <deal.II/lac/affine_constraints.h>
#include <deal.II/lac/trilinos_sparse_matrix.h>
#include <deal.II/lac/trilinos_parallel_block_vector.h>
#include <deal.II/lac/trilinos_precondition.h>
#include <deal.II/lac/solver_gmres.h>
#include <deal.II/lac/trilinos_solver.h>

#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/matrix_tools.h>

#include <fstream>
#include <iostream>
#include <memory>

using namespace dealii;

class FSIProblem
{
public:
  static constexpr unsigned int dim = 2;

  // Physical parameters
  struct Parameters
  {
    double viscosity = 1.0;
    double lambda = 1.0;
    double mu = 1.0;
  };

  // Boundary and Material IDs from .geo
  struct DomainLabels
  {
    static constexpr types::material_id fluid = 12;
    static constexpr types::material_id solid = 13;
    
    static constexpr types::boundary_id inflow = 14;
    static constexpr types::boundary_id bottom_clamped = 16;
  };

  class InletFunction : public Function<dim>
  {
  public:
    InletFunction() : Function<dim>(dim + 1 + dim) {} // [u, p, d]

    virtual void vector_value(const Point<dim> &p, Vector<double> &values) const override
    {
      values = 0.0;
      // Parabolic inlet on the y-component (matching typical FSI benchmarks)
      values[1] = -3.0 * std::sin(M_PI * (2.0 * p[0])); 
    }
  };

  FSIProblem(const std::string &mesh_file,
             const unsigned int velocity_degree,
             const unsigned int pressure_degree,
             const unsigned int displacement_degree);

  void setup();
  void assemble();
  void solve();
  void output();

private:
  // MPI context
  const unsigned int mpi_rank;
  const unsigned int mpi_size;
  ConditionalOStream pcout;

  // Discretization settings
  const std::string mesh_file_name;
  const unsigned int deg_u;
  const unsigned int deg_p;
  const unsigned int deg_d;

  Parameters physical_params;

  parallel::fullydistributed::Triangulation<dim> triangulation;
  DoFHandler<dim> dof_handler;

  // hp-FEM Collections
  hp::FECollection<dim> fe_collection;
  hp::QCollection<dim> quad_collection;
  hp::QCollection<dim - 1> face_quad_collection;
  hp::MappingCollection<dim> mapping_collection;

  AffineConstraints<double> system_constraints;

  // Linear algebra
  IndexSet locally_owned_dofs;
  IndexSet locally_relevant_dofs;

  TrilinosWrappers::SparseMatrix system_matrix;
  TrilinosWrappers::MPI::Vector system_rhs;
  TrilinosWrappers::MPI::Vector solution_owned;
  TrilinosWrappers::MPI::Vector solution_relevant;
};

#endif
