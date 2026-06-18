#include "FSIProblem.hpp"

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <iomanip>
#include <limits>

namespace
{
  // File-local aliases and helpers used by assembly, solves, and reporting.
  using MPIBlockMatrix = TrilinosWrappers::BlockSparseMatrix;
  using MPIBlockVector = TrilinosWrappers::MPI::BlockVector;
  using MPIVector = TrilinosWrappers::MPI::Vector;

  // Compact container for the diagnostics written after the solve.
  struct FieldStatistics
  {
    double l2_norm = 0.0;
    double min_value = std::numeric_limits<double>::infinity();
    double max_value = -std::numeric_limits<double>::infinity();
  };

  // Basic MPI reductions reused by the diagnostics code.
  double mpi_reduce_sum(const double local_value, const MPI_Comm communicator)
  {
    double global_value = 0.0;
    MPI_Allreduce(&local_value,
                  &global_value,
                  1,
                  MPI_DOUBLE,
                  MPI_SUM,
                  communicator);
    return global_value;
  }

  double mpi_reduce_min(const double local_value, const MPI_Comm communicator)
  {
    double global_value = 0.0;
    MPI_Allreduce(&local_value,
                  &global_value,
                  1,
                  MPI_DOUBLE,
                  MPI_MIN,
                  communicator);
    return global_value;
  }

  double mpi_reduce_max(const double local_value, const MPI_Comm communicator)
  {
    double global_value = 0.0;
    MPI_Allreduce(&local_value,
                  &global_value,
                  1,
                  MPI_DOUBLE,
                  MPI_MAX,
                  communicator);
    return global_value;
  }

  bool point_lexicographically_less(const Point<FSIProblem::dim> &lhs,
                                    const Point<FSIProblem::dim> &rhs)
  {
    // Keep rank-independent ordering when selecting or writing reference data.
    // Different processors (ranks) might compute the same point's coordinates but arrive at values that differ by a tiny fraction of a decimal.
    constexpr double tolerance = 1e-12;

    for (unsigned int d = 0; d < FSIProblem::dim; ++d)
    {
      if (lhs[d] < rhs[d] - tolerance)
        return true;

      if (lhs[d] > rhs[d] + tolerance)
        return false;
    }

    return false;
  }

  std::vector<unsigned int>
  match_face_quadrature_points(const FEFaceValues<FSIProblem::dim> &reference,
                               const FEFaceValues<FSIProblem::dim> &candidate)
  {
    const unsigned int n_q = reference.n_quadrature_points;

    AssertThrow(candidate.n_quadrature_points == n_q,
                ExcMessage("Interface faces use different quadrature sizes."));

    std::vector<unsigned int> permutation(n_q,
                                          std::numeric_limits<unsigned int>::max());
    std::vector<bool> used(n_q, false);

    constexpr double tolerance = 1e-12;
    for (unsigned int q = 0; q < n_q; ++q)
    {
      unsigned int best_index = std::numeric_limits<unsigned int>::max();
      double best_distance = std::numeric_limits<double>::infinity();

      for (unsigned int r = 0; r < n_q; ++r)
      {
        if (used[r])
          continue;

        const double distance =
            reference.quadrature_point(q).distance(candidate.quadrature_point(r));
        if (distance < best_distance)
        {
          best_distance = distance;
          best_index = r;
        }
      }

      AssertThrow(best_index != std::numeric_limits<unsigned int>::max() &&
                      best_distance < tolerance,
                  ExcMessage("Could not match fluid and solid interface "
                             "quadrature points."));

      permutation[q] = best_index;
      used[best_index] = true;
    }

    return permutation;
  }

  class StokesBlockOperator
  {
  public:
    explicit StokesBlockOperator(const MPIBlockMatrix &matrix_)
        : matrix(&matrix_)
    {
    }

    void vmult(MPIBlockVector &dst, const MPIBlockVector &src) const
    {
      // Apply the fluid [u, p] block without touching the solid unknowns.
      dst = 0.0;

      matrix->block(0, 0).vmult(dst.block(0), src.block(0));
      matrix->block(0, 1).vmult_add(dst.block(0), src.block(1));

      matrix->block(1, 0).vmult(dst.block(1), src.block(0));
      if (matrix->block(1, 1).n() != 0)
        matrix->block(1, 1).vmult_add(dst.block(1), src.block(1));
    }

  private:
    const MPIBlockMatrix *matrix;
  };

  class StokesBlockPreconditioner
  {
  public:
    StokesBlockPreconditioner(const MPIBlockMatrix &matrix_,
                              const TrilinosWrappers::SparseMatrix &pressure_mass_matrix_,
                              const std::vector<IndexSet> &block_owned_dofs_,
                              const MPI_Comm &mpi_communicator_,
                              const unsigned int degree_velocity)
        : matrix(&matrix_),
          pressure_mass_matrix(&pressure_mass_matrix_),
          tmp_velocity(block_owned_dofs_[0], mpi_communicator_)
    {
      TrilinosWrappers::PreconditionAMG::AdditionalData amg_data;
      amg_data.elliptic = true;
      amg_data.higher_order_elements = (degree_velocity > 1);
      amg_data.smoother_sweeps = 2;
      amg_data.aggregation_threshold = 1e-2;
      velocity_preconditioner.initialize(matrix->block(0, 0), amg_data);

      TrilinosWrappers::PreconditionJacobi::AdditionalData jacobi_data;
      pressure_preconditioner.initialize(*pressure_mass_matrix, jacobi_data);
    }

    void vmult(MPIBlockVector &dst, const MPIBlockVector &src) const
    {
      dst = 0.0;
      tmp_velocity = 0.0;

      // Approximate the Schur complement inverse with the pressure mass matrix.
      const double pressure_tolerance =
          std::max(1e-12, 1e-3 * src.block(1).l2_norm());
      SolverControl pressure_control(200, pressure_tolerance);
      SolverCG<MPIVector> pressure_solver(pressure_control);
      pressure_solver.solve(*pressure_mass_matrix,
                            dst.block(1),
                            src.block(1),
                            pressure_preconditioner);

      matrix->block(0, 1).vmult(tmp_velocity, dst.block(1));
      tmp_velocity *= -1.0;
      tmp_velocity += src.block(0);

      // Recover the velocity block once the pressure correction is available.
      const double velocity_tolerance =
          std::max(1e-12, 1e-3 * tmp_velocity.l2_norm());
      SolverControl velocity_control(200, velocity_tolerance);
      SolverCG<MPIVector> velocity_solver(velocity_control);
      velocity_solver.solve(matrix->block(0, 0),
                            dst.block(0),
                            tmp_velocity,
                            velocity_preconditioner);
    }

  private:
    const MPIBlockMatrix *matrix;
    const TrilinosWrappers::SparseMatrix *pressure_mass_matrix;

    TrilinosWrappers::PreconditionAMG velocity_preconditioner;
    TrilinosWrappers::PreconditionJacobi pressure_preconditioner;

    mutable MPIVector tmp_velocity;
  };
} // namespace.

void FSIProblem::setup()
{
  // Create the mesh.
  {
    pcout << "Initializing the mesh" << std::endl;

    // Read the mesh from file.
    GridIn<dim> grid_in;
    grid_in.attach_triangulation(mesh);

    std::ifstream mesh_file(mesh_file_name);
    AssertThrow(mesh_file.is_open(),
                ExcMessage("Cannot open mesh file '" + mesh_file_name +
                           "' from working directory '" +
                           std::filesystem::current_path().string() +
                           "'. Generate it with `cd ../mesh && make`."));

    unsigned int n_nodes = std::numeric_limits<unsigned int>::max();
    unsigned int n_elements = std::numeric_limits<unsigned int>::max();
    std::string token;
    while (mesh_file >> token)
    {
      if (token == "$Nodes")
        mesh_file >> n_nodes;
      else if (token == "$Elements")
        mesh_file >> n_elements;
    }

    AssertThrow(n_nodes != std::numeric_limits<unsigned int>::max(),
                ExcMessage("Mesh file '" + mesh_file_name +
                           "' does not contain a `$Nodes` section."));
    AssertThrow(n_elements != std::numeric_limits<unsigned int>::max(),
                ExcMessage("Mesh file '" + mesh_file_name +
                           "' does not contain an `$Elements` section."));
    AssertThrow(n_nodes > 0 && n_elements > 0,
                ExcMessage("Mesh file '" + mesh_file_name +
                           "' has zero nodes/elements. Regenerate it with "
                           "`cd ../mesh && make` (or `gmsh -2 mesh_fsi_geo.geo "
                           "-format msh2 -o mesh_fsi.msh`)."));

    mesh_file.clear();
    mesh_file.seekg(0);
    grid_in.read_msh(mesh_file);
    mesh.refine_global(2);

    pcout << "  Number of elements = " << mesh.n_active_cells() << std::endl;
    pcout << "  Number of MPI ranks = " << mpi_size << std::endl;
  }

  pcout << "-----------------------------------------------" << std::endl;

  // Initialize the finite element space.
  {
    pcout << "Initializing the finite element space" << std::endl;

    const FE_SimplexP<dim> fe_velocity(degree_velocity);
    const FE_SimplexP<dim> fe_pressure(degree_pressure);
    const auto rc = fe_velocity.reference_cell();
    const FE_Nothing<dim> fe_nothing(rc); // Used to fill the unused components of the block system.

    fe_stokes = std::make_unique<FESystem<dim>>(fe_velocity, dim,
                                                fe_pressure, 1,
                                                fe_nothing, dim);

    const FE_SimplexP<dim> fe_displacement(degree_displacement);
    fe_elasticity = std::make_unique<FESystem<dim>>(fe_nothing, dim,
                                                    FE_Nothing<dim>(rc), 1,
                                                    fe_displacement, dim);

    fe_collection.push_back(*fe_stokes);
    fe_collection.push_back(*fe_elasticity);

    mapping_collection.push_back(MappingFE<dim>(FE_SimplexP<dim>(1)));
    mapping_collection.push_back(MappingFE<dim>(FE_SimplexP<dim>(1)));

    quadrature_collection.push_back(QGaussSimplex<dim>(degree_velocity + 1));
    quadrature_collection.push_back(QGaussSimplex<dim>(degree_displacement + 1));

    face_quadrature_collection.push_back(QGaussSimplex<dim - 1>(degree_velocity + 1));
    face_quadrature_collection.push_back(QGaussSimplex<dim - 1>(degree_displacement + 1));

    pcout << "FE created" << std::endl;
    pcout << "fe_collection size = " << fe_collection.size() << std::endl;
    pcout << "  Velocity degree:           = " << fe_velocity.degree << std::endl;
    pcout << "  Pressure degree:           = " << fe_pressure.degree << std::endl;
    pcout << "  Displacement degree:       = " << fe_displacement.degree << std::endl;
    pcout << "  FE stokes components = " << fe_stokes->n_components() << std::endl;
    pcout << "  FE elasticity components = " << fe_elasticity->n_components()
          << std::endl;
  }

  pcout << "-----------------------------------------------" << std::endl;

  // Initialize the DoF handler.
  {
    pcout << "Initializing the DoF handler" << std::endl;

    dof_handler.reinit(mesh);

    for (const auto &cell : dof_handler.active_cell_iterators())
    {
      if (cell->material_id() == MaterialIds::Fluid)
        cell->set_active_fe_index(0);
      else if (cell->material_id() == MaterialIds::Solid)
        cell->set_active_fe_index(1);
      else
        AssertThrow(false,
                    ExcMessage("Unexpected material id " +
                               std::to_string(cell->material_id()) +
                               ". Expected Fluid=12 or Solid=13."));
    }

    dof_handler.distribute_dofs(fe_collection);

    // Reorder unknowns into the global block layout [u, p, d].
    std::vector<unsigned int> block_component(2 * dim + 1, 0);
    block_component[dim] = 1;
    for (unsigned int c = dim + 1; c < 2 * dim + 1; ++c)
      block_component[c] = 2;

    DoFRenumbering::component_wise(dof_handler, block_component);

    locally_owned_dofs = dof_handler.locally_owned_dofs();
    DoFTools::extract_locally_relevant_dofs(dof_handler, locally_relevant_dofs);

    const std::vector<types::global_dof_index> dofs_per_block =
        DoFTools::count_dofs_per_fe_block(dof_handler, block_component);

    block_offsets.resize(4);
    block_offsets[0] = 0;
    block_offsets[1] = dofs_per_block[0];
    block_offsets[2] = dofs_per_block[0] + dofs_per_block[1];
    block_offsets[3] = dof_handler.n_dofs();

    block_owned_dofs.resize(3);
    block_relevant_dofs.resize(3);

    for (unsigned int b = 0; b < 3; ++b)
    {
      block_owned_dofs[b] =
          locally_owned_dofs.get_view(block_offsets[b], block_offsets[b + 1]);
      block_relevant_dofs[b] =
          locally_relevant_dofs.get_view(block_offsets[b], block_offsets[b + 1]);
    }

    pcout << "DoF: done distribute" << std::endl;
    pcout << "Total DoFs = " << dof_handler.n_dofs() << std::endl;
    pcout << "  Block DoFs: velocity = " << block_offsets[1]
          << " pressure = " << (block_offsets[2] - block_offsets[1])
          << " displacement = " << (block_offsets[3] - block_offsets[2])
          << std::endl;
    pcout << "  Locally owned on rank " << mpi_rank
          << ": velocity = " << block_owned_dofs[0].n_elements()
          << ", pressure = " << block_owned_dofs[1].n_elements()
          << ", displacement = " << block_owned_dofs[2].n_elements() << std::endl;
  }

  pcout << "-----------------------------------------------" << std::endl;

  // Initialize the linear system.
  {
    pcout << "Initializing the linear system" << std::endl;
    pcout << "  Initializing the sparsity pattern" << std::endl;

    Table<2, DoFTools::Coupling> cell_coupling(2 * dim + 1, 2 * dim + 1);
    Table<2, DoFTools::Coupling> face_coupling(2 * dim + 1, 2 * dim + 1);

    for (unsigned int c = 0; c < 2 * dim + 1; ++c)
      for (unsigned int d = 0; d < 2 * dim + 1; ++d)
      {
        cell_coupling[c][d] = DoFTools::none;
        face_coupling[c][d] = DoFTools::none;
      }

    for (unsigned int c = 0; c < dim; ++c)
      for (unsigned int d = 0; d < dim; ++d)
        cell_coupling[c][d] = DoFTools::always;

    for (unsigned int c = 0; c < dim; ++c)
    {
      cell_coupling[c][dim] = DoFTools::always;
      cell_coupling[dim][c] = DoFTools::always;
    }

    for (unsigned int c = dim + 1; c < 2 * dim + 1; ++c)
      for (unsigned int d = dim + 1; d < 2 * dim + 1; ++d)
        cell_coupling[c][d] = DoFTools::always;

    for (unsigned int c = dim + 1; c < 2 * dim + 1; ++c)
      for (unsigned int d = 0; d <= dim; ++d)
        face_coupling[c][d] = DoFTools::always;

    constraints.clear(); // velocity is null (0) on the interface with constraints.
    constraints.reinit(locally_relevant_dofs);

    DoFTools::make_hanging_node_constraints(dof_handler, constraints);

    // Constrain the fluid velocity DoFs that lie on the fluid-solid interface.
    std::vector<types::global_dof_index> local_face_dof_indices(fe_stokes->n_dofs_per_face());

    for (const auto &cell : dof_handler.active_cell_iterators())
    {
      if (cell->active_fe_index() != 0)
        continue;

      for (const auto face_no : cell->face_indices())
      {
        if (cell->face(face_no)->at_boundary())
          continue;

        bool face_is_on_interface = false;
        const auto neighbor = cell->neighbor(face_no);

        if (!neighbor->has_children() && neighbor->active_fe_index() == 1)
          face_is_on_interface = true;
        else if (neighbor->has_children())
          for (unsigned int sf = 0; sf < cell->face(face_no)->n_children(); ++sf)
            if (cell->neighbor_child_on_subface(face_no, sf)->active_fe_index() == 1)
            {
              face_is_on_interface = true;
              break;
            }

        if (!face_is_on_interface)
          continue;

        cell->face(face_no)->get_dof_indices(local_face_dof_indices, 0);

        for (unsigned int i = 0; i < local_face_dof_indices.size(); ++i)
        {
          const types::global_dof_index dof_index = local_face_dof_indices[i];

          if (fe_stokes->face_system_to_component_index(i).first < dim &&
              locally_relevant_dofs.is_element(dof_index))
          {
            if (!constraints.is_constrained(dof_index))
              constraints.add_line(dof_index);

            constraints.set_inhomogeneity(dof_index, 0.0);
          }
        }
      }
    }

    ComponentMask mask_velocity(2 * dim + 1, false);
    for (unsigned int c = 0; c < dim; ++c)
      mask_velocity.set(c, true);

    ComponentMask mask_displacement(2 * dim + 1, false);
    for (unsigned int c = dim + 1; c < 2 * dim + 1; ++c)
      mask_displacement.set(c, true);

    Functions::ZeroFunction<dim> zero_function(2 * dim + 1);

    VectorTools::interpolate_boundary_values(dof_handler,
                                             BoundaryIds::Inflow,
                                             inlet_velocity,
                                             constraints,
                                             mask_velocity);
    // VectorTools::interpolate_boundary_values(dof_handler,
    //                                         BoundaryIds::LeftWallFluid,
    //                                         zero_function,
    //                                         constraints,
    //                                         mask_velocity);
    // VectorTools::interpolate_boundary_values(dof_handler,
    //                                         BoundaryIds::RightWallFluid,
    //                                         zero_function,
    //                                         constraints,
    //                                         mask_velocity);
    VectorTools::interpolate_boundary_values(dof_handler,
                                             BoundaryIds::BottomWall,
                                             zero_function,
                                             constraints,
                                             mask_displacement);

    // Pin one pressure DoF in a deterministic location to remove the null space.
    types::global_dof_index pressure_reference =
        numbers::invalid_dof_index;
    Point<dim> pressure_reference_center;
    for (unsigned int d = 0; d < dim; ++d)
      pressure_reference_center[d] = std::numeric_limits<double>::infinity();

    std::vector<types::global_dof_index> candidate_dof_indices;
    for (const auto &cell : dof_handler.active_cell_iterators())
    {
      if (cell->active_fe_index() != 0)
        continue;

      const Point<dim> cell_center = cell->center();
      if (pressure_reference != numbers::invalid_dof_index &&
          !point_lexicographically_less(cell_center, pressure_reference_center))
        continue;

      candidate_dof_indices.resize(cell->get_fe().n_dofs_per_cell());
      cell->get_dof_indices(candidate_dof_indices);

      for (unsigned int i = 0; i < candidate_dof_indices.size(); ++i)
        if (cell->get_fe().system_to_component_index(i).first == dim)
        {
          pressure_reference = candidate_dof_indices[i];
          pressure_reference_center = cell_center;
          break;
        }
    }

    AssertThrow(pressure_reference != numbers::invalid_dof_index,
                ExcMessage("Could not find a deterministic pressure reference DoF."));

    if (locally_relevant_dofs.is_element(pressure_reference))
    {
      if (!constraints.is_constrained(pressure_reference))
        constraints.add_line(pressure_reference);
      constraints.set_inhomogeneity(pressure_reference, 0.0);
    }
    pcout << "  Pressure reference DoF = " << pressure_reference
          << " at (" << pressure_reference_center[0] << ", "
          << pressure_reference_center[1] << ")" << std::endl;

    constraints.close();

    sparsity.reinit(block_owned_dofs,
                    block_owned_dofs,
                    block_relevant_dofs,
                    mpi_communicator);
    DoFTools::make_flux_sparsity_pattern(dof_handler,
                                         sparsity,
                                         cell_coupling,
                                         face_coupling);
    for (const auto index : locally_relevant_dofs)
      if (constraints.is_constrained(index))
        sparsity.add(index, index);
    sparsity.compress();

    DynamicSparsityPattern pressure_mass_dsp(block_relevant_dofs[1]);
    // Build the pressure mass matrix separately for the Stokes preconditioner.
    for (const auto &cell : dof_handler.active_cell_iterators())
    {
      if (!cell->is_locally_owned() || cell->active_fe_index() != 0)
        continue;

      std::vector<types::global_dof_index> local_dof_indices(cell->get_fe().n_dofs_per_cell());
      cell->get_dof_indices(local_dof_indices);

      for (unsigned int i = 0; i < local_dof_indices.size(); ++i)
      {
        if (cell->get_fe().system_to_component_index(i).first != dim)
          continue;

        const auto pressure_i = local_dof_indices[i] - block_offsets[1];
        for (unsigned int j = 0; j < local_dof_indices.size(); ++j)
        {
          if (cell->get_fe().system_to_component_index(j).first != dim)
            continue;

          const auto pressure_j = local_dof_indices[j] - block_offsets[1];
          pressure_mass_dsp.add(pressure_i, pressure_j);
        }
      }
    }

    SparsityTools::distribute_sparsity_pattern(pressure_mass_dsp,
                                               block_owned_dofs[1],
                                               mpi_communicator,
                                               block_relevant_dofs[1]);

    TrilinosWrappers::SparsityPattern pressure_mass_sparsity;
    pressure_mass_sparsity.reinit(block_owned_dofs[1],
                                  block_owned_dofs[1],
                                  pressure_mass_dsp,
                                  mpi_communicator);
    pressure_mass_sparsity.compress();

    pcout << "  Initializing the matrices" << std::endl;
    system_matrix.reinit(sparsity);
    pressure_mass_matrix.reinit(pressure_mass_sparsity);

    pcout << "  Initializing the system right-hand side" << std::endl;
    system_rhs.reinit(block_owned_dofs, mpi_communicator);
    pcout << "  Initializing the solution vector" << std::endl;
    solution_owned.reinit(block_owned_dofs, mpi_communicator);
    solution.reinit(block_owned_dofs, block_relevant_dofs, mpi_communicator);
  }
}

void FSIProblem::assemble()
{
  pcout << "===============================================" << std::endl;
  pcout << "Assembling the system" << std::endl;

  hp::FEValues<dim> hp_fe_values(mapping_collection,
                                 fe_collection,
                                 quadrature_collection,
                                 update_values | update_gradients |
                                     update_quadrature_points | update_JxW_values);
  const UpdateFlags face_update_flags =
      update_values | update_gradients | update_quadrature_points |
      update_normal_vectors | update_JxW_values;
  hp::FEFaceValues<dim> hp_fe_face_values_solid(mapping_collection,
                                                fe_collection,
                                                face_quadrature_collection,
                                                face_update_flags);
  hp::FEFaceValues<dim> hp_fe_face_values_fluid(mapping_collection,
                                                fe_collection,
                                                face_quadrature_collection,
                                                face_update_flags);

  FEValuesExtractors::Vector velocity(0);
  FEValuesExtractors::Scalar pressure(dim);
  FEValuesExtractors::Vector displacement(dim + 1);

  system_matrix = 0.0;
  pressure_mass_matrix = 0.0;
  system_rhs = 0.0;

  std::vector<types::global_dof_index> dof_indices;

  for (const auto &cell : dof_handler.active_cell_iterators())
  {
    if (!cell->is_locally_owned())
      continue;

    const unsigned int fe_index = cell->active_fe_index();
    const unsigned int dofs_per_cell = cell->get_fe().n_dofs_per_cell();

    FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell);
    FullMatrix<double> cell_pressure_mass(dofs_per_cell, dofs_per_cell);
    Vector<double> cell_rhs(dofs_per_cell);

    cell_matrix = 0.0;
    cell_pressure_mass = 0.0;
    cell_rhs = 0.0;

    dof_indices.resize(dofs_per_cell);

    hp_fe_values.reinit(cell);
    const FEValues<dim> &fe_values = hp_fe_values.get_present_fe_values();

    if (fe_index == 0)
    {
      // Fluid cells contribute the Stokes block and the pressure mass matrix.
      for (unsigned int q = 0; q < fe_values.n_quadrature_points; ++q)
        for (unsigned int i = 0; i < dofs_per_cell; ++i)
          for (unsigned int j = 0; j < dofs_per_cell; ++j)
          {
            cell_matrix(i, j) +=
                nu *
                scalar_product(fe_values[velocity].gradient(i, q),
                               fe_values[velocity].gradient(j, q)) *
                fe_values.JxW(q);

            cell_matrix(i, j) -= fe_values[velocity].divergence(i, q) *
                                 fe_values[pressure].value(j, q) *
                                 fe_values.JxW(q);

            cell_matrix(i, j) -= fe_values[velocity].divergence(j, q) *
                                 fe_values[pressure].value(i, q) *
                                 fe_values.JxW(q);

            cell_pressure_mass(i, j) +=
                fe_values[pressure].value(i, q) *
                fe_values[pressure].value(j, q) / nu * fe_values.JxW(q);
          }
    }
    else if (fe_index == 1)
    {
      // Solid cells contribute the linear elasticity block.
      for (unsigned int q = 0; q < fe_values.n_quadrature_points; ++q)
        for (unsigned int i = 0; i < dofs_per_cell; ++i)
          for (unsigned int j = 0; j < dofs_per_cell; ++j)
            cell_matrix(i, j) +=
                (2.0 * mu *
                     scalar_product(fe_values[displacement].symmetric_gradient(j, q),
                                    fe_values[displacement].symmetric_gradient(i, q)) +
                 lambda * fe_values[displacement].divergence(j, q) *
                     fe_values[displacement].divergence(i, q)) *
                fe_values.JxW(q);
    }

    cell->get_dof_indices(dof_indices);
    constraints.distribute_local_to_global(cell_matrix,
                                           cell_rhs,
                                           dof_indices,
                                           system_matrix,
                                           system_rhs);

    if (fe_index == 0)
      for (unsigned int i = 0; i < dofs_per_cell; ++i)
      {
        if (cell->get_fe().system_to_component_index(i).first != dim)
          continue;

        const auto pressure_i = dof_indices[i] - block_offsets[1];
        for (unsigned int j = 0; j < dofs_per_cell; ++j)
        {
          if (cell->get_fe().system_to_component_index(j).first != dim)
            continue;

          const auto pressure_j = dof_indices[j] - block_offsets[1];
          pressure_mass_matrix.add(pressure_i,
                                   pressure_j,
                                   cell_pressure_mass(i, j));
        }
      }

    if (fe_index != 1)
      continue;

    // Add the fluid traction acting on the solid side of each interface face.
    for (unsigned int f = 0; f < cell->n_faces(); ++f)
    {
      if (cell->face(f)->at_boundary())
        continue;

      const auto neighbor = cell->neighbor(f);
      if (neighbor->active_fe_index() != 0)
        continue;

      hp_fe_face_values_solid.reinit(cell, f);
      const FEFaceValues<dim> &fe_face_solid =
          hp_fe_face_values_solid.get_present_fe_values();

      const unsigned int f_fluid = cell->neighbor_of_neighbor(f);
      hp_fe_face_values_fluid.reinit(neighbor, f_fluid);
      const FEFaceValues<dim> &fe_face_fluid =
          hp_fe_face_values_fluid.get_present_fe_values();
      const std::vector<unsigned int> fluid_quadrature_index =
          match_face_quadrature_points(fe_face_solid, fe_face_fluid);

      const unsigned int dofs_solid = cell->get_fe().dofs_per_cell;
      const unsigned int dofs_fluid = neighbor->get_fe().dofs_per_cell;

      FullMatrix<double> local_interface_matrix(dofs_solid, dofs_fluid);
      local_interface_matrix = 0.0;

      for (unsigned int q = 0; q < fe_face_solid.n_quadrature_points; ++q)
      {
        const Tensor<1, dim> n_solid = fe_face_solid.normal_vector(q);

        for (unsigned int i = 0; i < dofs_solid; ++i)
        {
          const Tensor<1, dim> phi_d_i =
              fe_face_solid[displacement].value(i, q);

          for (unsigned int j = 0; j < dofs_fluid; ++j)
          {
            const unsigned int q_fluid = fluid_quadrature_index[q];
            const Tensor<1, dim> traction =
                -(nu * (fe_face_fluid[velocity].gradient(j, q_fluid) * n_solid) -
                  fe_face_fluid[pressure].value(j, q_fluid) * n_solid);

            local_interface_matrix(i, j) +=
                (traction * phi_d_i) * fe_face_solid.JxW(q);
          }
        }
      }

      std::vector<types::global_dof_index> dofs_solid_indices(dofs_solid);
      std::vector<types::global_dof_index> dofs_fluid_indices(dofs_fluid);
      cell->get_dof_indices(dofs_solid_indices);
      neighbor->get_dof_indices(dofs_fluid_indices);

      constraints.distribute_local_to_global(local_interface_matrix,
                                             dofs_solid_indices,
                                             dofs_fluid_indices,
                                             system_matrix);
    }
  }

  system_matrix.compress(VectorOperation::add);
  pressure_mass_matrix.compress(VectorOperation::add);
  system_rhs.compress(VectorOperation::add);
}

void FSIProblem::solve()
{
  pcout << "===============================================" << std::endl;
  pcout << "Solving the linear system" << std::endl;

  solution_owned = 0.0;

  // Solve the fluid [u, p] block first with a block preconditioned Krylov method.
  std::vector<IndexSet> stokes_owned(2), stokes_relevant(2);
  stokes_owned[0] = block_owned_dofs[0];
  stokes_owned[1] = block_owned_dofs[1];
  stokes_relevant[0] = block_relevant_dofs[0];
  stokes_relevant[1] = block_relevant_dofs[1];

  MPIBlockVector stokes_rhs;
  MPIBlockVector stokes_solution;
  stokes_rhs.reinit(stokes_owned, mpi_communicator);
  stokes_solution.reinit(stokes_owned, mpi_communicator);
  stokes_rhs.block(0) = system_rhs.block(0);
  stokes_rhs.block(1) = system_rhs.block(1);

  StokesBlockOperator stokes_matrix(system_matrix);
  StokesBlockPreconditioner stokes_preconditioner(system_matrix,
                                                  pressure_mass_matrix,
                                                  block_owned_dofs,
                                                  mpi_communicator,
                                                  degree_velocity);

  const double fluid_rhs_norm = stokes_rhs.l2_norm();
  SolverControl stokes_control(400, std::max(1e-12, 1e-8 * fluid_rhs_norm));
  SolverFGMRES<MPIBlockVector> stokes_solver(stokes_control);
  stokes_solver.solve(stokes_matrix,
                      stokes_solution,
                      stokes_rhs,
                      stokes_preconditioner);

  solution_owned.block(0) = stokes_solution.block(0);
  solution_owned.block(1) = stokes_solution.block(1);

  // Move the fluid-induced interface terms to the solid right-hand side.
  MPIVector solid_rhs(block_owned_dofs[2], mpi_communicator);
  MPIVector solid_tmp(block_owned_dofs[2], mpi_communicator);
  solid_rhs = system_rhs.block(2);

  system_matrix.block(2, 0).vmult(solid_tmp, solution_owned.block(0));
  solid_rhs -= solid_tmp;
  system_matrix.block(2, 1).vmult(solid_tmp, solution_owned.block(1));
  solid_rhs -= solid_tmp;

  TrilinosWrappers::PreconditionAMG::AdditionalData solid_amg_data;
  solid_amg_data.elliptic = true;
  solid_amg_data.higher_order_elements = (degree_displacement > 1);
  solid_amg_data.smoother_sweeps = 8;
  solid_amg_data.aggregation_threshold = 5e-2;

  TrilinosWrappers::PreconditionAMG solid_preconditioner;
  solid_preconditioner.initialize(system_matrix.block(2, 2), solid_amg_data);

  const double solid_rhs_norm = solid_rhs.l2_norm();
  SolverControl solid_control(400, std::max(1e-12, 1e-8 * solid_rhs_norm));
  SolverCG<MPIVector> solid_solver(solid_control);
  solid_solver.solve(system_matrix.block(2, 2),
                     solution_owned.block(2),
                     solid_rhs,
                     solid_preconditioner);

  // Enforce constrained values and refresh the ghosted view used for output.
  constraints.distribute(solution_owned);

  solution = solution_owned;
  solution.update_ghost_values();

  last_stokes_iterations = stokes_control.last_step();
  last_solid_iterations = solid_control.last_step();

  pcout << "  Stokes iterations = " << last_stokes_iterations << std::endl;
  pcout << "  Solid iterations  = " << last_solid_iterations << std::endl;
}

void FSIProblem::write_diagnostics_report(const std::string &output_file_name) const
{
  // Compute global norms and extrema from cell-wise quadrature samples.
  hp::FEValues<dim> hp_fe_values(mapping_collection,
                                 fe_collection,
                                 quadrature_collection,
                                 update_values | update_JxW_values);

  FEValuesExtractors::Vector velocity(0);
  FEValuesExtractors::Scalar pressure(dim);
  FEValuesExtractors::Vector displacement(dim + 1);

  FieldStatistics velocity_statistics;
  FieldStatistics pressure_statistics;
  FieldStatistics displacement_statistics;

  double local_velocity_l2_sq = 0.0;
  double local_pressure_l2_sq = 0.0;
  double local_displacement_l2_sq = 0.0;

  for (const auto &cell : dof_handler.active_cell_iterators())
  {
    if (!cell->is_locally_owned())
      continue;

    hp_fe_values.reinit(cell);
    const FEValues<dim> &fe_values = hp_fe_values.get_present_fe_values();
    const unsigned int n_q = fe_values.n_quadrature_points;

    if (cell->active_fe_index() == 0)
    {
      std::vector<Tensor<1, dim>> velocity_values(n_q);
      std::vector<double> pressure_values(n_q);

      fe_values[velocity].get_function_values(solution, velocity_values);
      fe_values[pressure].get_function_values(solution, pressure_values);

      for (unsigned int q = 0; q < n_q; ++q)
      {
        const double velocity_magnitude = velocity_values[q].norm();
        local_velocity_l2_sq +=
            velocity_values[q].norm_square() * fe_values.JxW(q);
        local_pressure_l2_sq +=
            pressure_values[q] * pressure_values[q] * fe_values.JxW(q);

        velocity_statistics.min_value =
            std::min(velocity_statistics.min_value, velocity_magnitude);
        velocity_statistics.max_value =
            std::max(velocity_statistics.max_value, velocity_magnitude);

        pressure_statistics.min_value =
            std::min(pressure_statistics.min_value, pressure_values[q]);
        pressure_statistics.max_value =
            std::max(pressure_statistics.max_value, pressure_values[q]);
      }
    }
    else if (cell->active_fe_index() == 1)
    {
      std::vector<Tensor<1, dim>> displacement_values(n_q);
      fe_values[displacement].get_function_values(solution, displacement_values);

      for (unsigned int q = 0; q < n_q; ++q)
      {
        const double displacement_magnitude =
            displacement_values[q].norm();
        local_displacement_l2_sq +=
            displacement_values[q].norm_square() * fe_values.JxW(q);

        displacement_statistics.min_value =
            std::min(displacement_statistics.min_value,
                     displacement_magnitude);
        displacement_statistics.max_value =
            std::max(displacement_statistics.max_value,
                     displacement_magnitude);
      }
    }
  }

  velocity_statistics.l2_norm =
      std::sqrt(mpi_reduce_sum(local_velocity_l2_sq, mpi_communicator));
  pressure_statistics.l2_norm =
      std::sqrt(mpi_reduce_sum(local_pressure_l2_sq, mpi_communicator));
  displacement_statistics.l2_norm =
      std::sqrt(mpi_reduce_sum(local_displacement_l2_sq, mpi_communicator));

  velocity_statistics.min_value =
      mpi_reduce_min(velocity_statistics.min_value, mpi_communicator);
  velocity_statistics.max_value =
      mpi_reduce_max(velocity_statistics.max_value, mpi_communicator);
  pressure_statistics.min_value =
      mpi_reduce_min(pressure_statistics.min_value, mpi_communicator);
  pressure_statistics.max_value =
      mpi_reduce_max(pressure_statistics.max_value, mpi_communicator);
  displacement_statistics.min_value =
      mpi_reduce_min(displacement_statistics.min_value, mpi_communicator);
  displacement_statistics.max_value =
      mpi_reduce_max(displacement_statistics.max_value, mpi_communicator);

  pcout << "  Velocity |u|_L2       = " << velocity_statistics.l2_norm
        << "  min/max |u| = [" << velocity_statistics.min_value << ", "
        << velocity_statistics.max_value << "]" << std::endl;
  pcout << "  Pressure ||p||_L2     = " << pressure_statistics.l2_norm
        << "  min/max p = [" << pressure_statistics.min_value << ", "
        << pressure_statistics.max_value << "]" << std::endl;
  pcout << "  Displacement |d|_L2   = " << displacement_statistics.l2_norm
        << "  min/max |d| = [" << displacement_statistics.min_value << ", "
        << displacement_statistics.max_value << "]" << std::endl;

  if (mpi_rank != 0)
    return;

  // Emit a compact rank-0 report that can be consumed by scripts.
  std::ofstream report(output_file_name + "_np" + std::to_string(mpi_size) +
                       ".report");
  AssertThrow(report.is_open(),
              ExcMessage("Could not open diagnostics report file for writing."));

  report << std::setprecision(17);
  report << "mpi_ranks=" << mpi_size << '\n';
  report << "n_dofs=" << dof_handler.n_dofs() << '\n';
  report << "stokes_iterations=" << last_stokes_iterations << '\n';
  report << "solid_iterations=" << last_solid_iterations << '\n';
  report << "velocity_l2=" << velocity_statistics.l2_norm << '\n';
  report << "velocity_magnitude_min=" << velocity_statistics.min_value << '\n';
  report << "velocity_magnitude_max=" << velocity_statistics.max_value << '\n';
  report << "pressure_l2=" << pressure_statistics.l2_norm << '\n';
  report << "pressure_min=" << pressure_statistics.min_value << '\n';
  report << "pressure_max=" << pressure_statistics.max_value << '\n';
  report << "displacement_l2=" << displacement_statistics.l2_norm << '\n';
  report << "displacement_magnitude_min="
         << displacement_statistics.min_value << '\n';
  report << "displacement_magnitude_max="
         << displacement_statistics.max_value << '\n';
  report << "pvtu_file=./" << output_file_name << "_0000.pvtu\n";
  report << "solution_snapshot_file=./" << output_file_name << "_np" << mpi_size
         << ".solution.csv\n";
}

void FSIProblem::write_solution_snapshot(const std::string &output_file_name) const
{
  // Export one deterministic, cell-averaged sample per active cell.
  hp::FEValues<dim> hp_fe_values(mapping_collection,
                                 fe_collection,
                                 quadrature_collection,
                                 update_values | update_JxW_values);

  FEValuesExtractors::Vector velocity(0);
  FEValuesExtractors::Scalar pressure(dim);
  FEValuesExtractors::Vector displacement(dim + 1);

  constexpr unsigned int row_size = 3 * dim + 3;

  std::vector<double> local_rows;
  local_rows.reserve(mesh.n_active_cells() * row_size);

  for (const auto &cell : dof_handler.active_cell_iterators())
  {
    if (!cell->is_locally_owned())
      continue;

    hp_fe_values.reinit(cell);
    const FEValues<dim> &fe_values = hp_fe_values.get_present_fe_values();
    const unsigned int n_q = fe_values.n_quadrature_points;
    const double cell_measure = cell->measure();

    Tensor<1, dim> velocity_average;
    double pressure_average = 0.0;
    Tensor<1, dim> displacement_average;

    if (cell->active_fe_index() == 0)
    {
      std::vector<Tensor<1, dim>> velocity_values(n_q);
      std::vector<double> pressure_values(n_q);

      fe_values[velocity].get_function_values(solution, velocity_values);
      fe_values[pressure].get_function_values(solution, pressure_values);

      for (unsigned int q = 0; q < n_q; ++q)
      {
        velocity_average += velocity_values[q] * fe_values.JxW(q);
        pressure_average += pressure_values[q] * fe_values.JxW(q);
      }

      velocity_average /= cell_measure;
      pressure_average /= cell_measure;
    }
    else if (cell->active_fe_index() == 1)
    {
      std::vector<Tensor<1, dim>> displacement_values(n_q);
      fe_values[displacement].get_function_values(solution, displacement_values);

      for (unsigned int q = 0; q < n_q; ++q)
        displacement_average += displacement_values[q] * fe_values.JxW(q);

      displacement_average /= cell_measure;
    }

    local_rows.push_back(static_cast<double>(cell->material_id()));
    for (unsigned int d = 0; d < dim; ++d)
      local_rows.push_back(cell->center()[d]);
    local_rows.push_back(cell_measure);
    for (unsigned int d = 0; d < dim; ++d)
      local_rows.push_back(velocity_average[d]);
    local_rows.push_back(pressure_average);
    for (unsigned int d = 0; d < dim; ++d)
      local_rows.push_back(displacement_average[d]);
  }

  const int local_count = static_cast<int>(local_rows.size());

  std::vector<int> counts;
  if (mpi_rank == 0)
    counts.resize(mpi_size);

  MPI_Gather(&local_count,
             1,
             MPI_INT,
             mpi_rank == 0 ? counts.data() : nullptr,
             1,
             MPI_INT,
             0,
             mpi_communicator);

  std::vector<int> displacements;
  int global_count = 0;
  if (mpi_rank == 0)
  {
    displacements.resize(mpi_size, 0);
    for (unsigned int rank = 0; rank < mpi_size; ++rank)
    {
      displacements[rank] = global_count;
      global_count += counts[rank];
    }
  }

  std::vector<double> gathered_rows;
  if (mpi_rank == 0)
    gathered_rows.resize(global_count);

  MPI_Gatherv(local_rows.data(),
              local_count,
              MPI_DOUBLE,
              mpi_rank == 0 ? gathered_rows.data() : nullptr,
              mpi_rank == 0 ? counts.data() : nullptr,
              mpi_rank == 0 ? displacements.data() : nullptr,
              MPI_DOUBLE,
              0,
              mpi_communicator);

  if (mpi_rank != 0)
    return;

  AssertThrow(global_count % static_cast<int>(row_size) == 0,
              ExcMessage("Gathered cell-sample snapshot has an invalid size."));

  // Sort by material and cell center so snapshots stay comparable across runs.
  std::vector<unsigned int> row_order(global_count / row_size);
  for (unsigned int i = 0; i < row_order.size(); ++i)
    row_order[i] = i;

  std::sort(row_order.begin(),
            row_order.end(),
            [&](const unsigned int lhs, const unsigned int rhs)
            {
              const unsigned int lhs_offset = lhs * row_size;
              const unsigned int rhs_offset = rhs * row_size;

              const auto lhs_material =
                  static_cast<unsigned int>(gathered_rows[lhs_offset]);
              const auto rhs_material =
                  static_cast<unsigned int>(gathered_rows[rhs_offset]);

              if (lhs_material != rhs_material)
                return lhs_material < rhs_material;

              Point<dim> lhs_center;
              Point<dim> rhs_center;
              for (unsigned int d = 0; d < dim; ++d)
              {
                lhs_center[d] = gathered_rows[lhs_offset + 1 + d];
                rhs_center[d] = gathered_rows[rhs_offset + 1 + d];
              }

              return point_lexicographically_less(lhs_center, rhs_center);
            });

  std::ofstream snapshot(output_file_name + "_np" + std::to_string(mpi_size) +
                         ".solution.csv");
  AssertThrow(snapshot.is_open(),
              ExcMessage("Could not open solution snapshot file for writing."));

  snapshot << std::setprecision(17);
  snapshot << "material,center_x,center_y,measure,"
              "velocity_x,velocity_y,pressure,"
              "displacement_x,displacement_y\n";

  for (const unsigned int row : row_order)
  {
    const unsigned int offset = row * row_size;
    snapshot << static_cast<unsigned int>(gathered_rows[offset]);
    for (unsigned int i = 1; i < row_size; ++i)
      snapshot << ',' << gathered_rows[offset + i];
    snapshot << '\n';
  }
}

void FSIProblem::output()
{
  pcout << "===============================================" << std::endl;

  DataOut<dim> data_out;
  data_out.attach_dof_handler(dof_handler);

  std::vector<std::string> names;
  std::vector<DataComponentInterpretation::DataComponentInterpretation> interpretation;

  for (unsigned int d = 0; d < dim; ++d)
  {
    names.push_back("velocity");
    interpretation.push_back(DataComponentInterpretation::component_is_part_of_vector);
  }

  names.push_back("pressure");
  interpretation.push_back(DataComponentInterpretation::component_is_scalar);

  for (unsigned int d = 0; d < dim; ++d)
  {
    names.push_back("displacement");
    interpretation.push_back(DataComponentInterpretation::component_is_part_of_vector);
  }

  data_out.add_data_vector(dof_handler, solution, names, interpretation);

  std::vector<unsigned int> partition_int(mesh.n_active_cells());
  GridTools::get_subdomain_association(mesh, partition_int);
  const dealii::Vector<double> partitioning(partition_int.begin(), partition_int.end());
  data_out.add_data_vector(partitioning, "partitioning");

  data_out.build_patches();

  const std::string output_file_name = "output-FSIProblem";
  // Write visualization data first, then the text-based diagnostics products.
  data_out.write_vtu_with_pvtu_record("./",
                                      output_file_name,
                                      0,
                                      mpi_communicator,
                                      4);
  write_diagnostics_report(output_file_name);
  write_solution_snapshot(output_file_name);

  pcout << "Output written to ./" << output_file_name << "_0000.pvtu"
        << std::endl;
  pcout << "Diagnostics written to ./" << output_file_name << "_np" << mpi_size
        << ".report" << std::endl;
  pcout << "Solution snapshot written to ./" << output_file_name << "_np"
        << mpi_size << ".solution.csv" << std::endl;
  pcout << "===============================================" << std::endl;
}
