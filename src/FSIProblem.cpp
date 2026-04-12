#include "FSIProblem.hpp"

FSIProblem::FSIProblem(const std::string &mesh_file,
                       const unsigned int velocity_degree,
                       const unsigned int pressure_degree,
                       const unsigned int displacement_degree)
  : mpi_rank(Utilities::MPI::this_mpi_process(MPI_COMM_WORLD))
  , mpi_size(Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD))
  , pcout(std::cout, mpi_rank == 0)
  , mesh_file_name(mesh_file)
  , deg_u(velocity_degree)
  , deg_p(pressure_degree)
  , deg_d(displacement_degree)
  , triangulation(MPI_COMM_WORLD)
{}

void FSIProblem::setup()
{
  pcout << "Building Mesh and Finite Element Collections" << std::endl;

  {
    Triangulation<dim> serial_triangulation;
    GridIn<dim> grid_in;
    grid_in.attach_triangulation(serial_triangulation);
    std::ifstream input_mesh(mesh_file_name);
    grid_in.read_msh(input_mesh);

    GridTools::partition_triangulation(mpi_size, serial_triangulation);
    const auto description = 
      TriangulationDescription::Utilities::create_description_from_triangulation(
        serial_triangulation, MPI_COMM_WORLD);
    triangulation.create_triangulation(description);
    pcout << "  Global active cells: " << triangulation.n_global_active_cells() << std::endl;
  }

  // 1. FE Collection setup
  {
    const FE_SimplexP<dim> fe_u(deg_u);
    const FE_SimplexP<dim> fe_p(deg_p);
    const FE_SimplexP<dim> fe_d(deg_d);
    
    // Components: [u, p, d]
    // Fluid system: u (dim), p (1), d (dim, deactivated)
    fe_collection.push_back(
      FESystem<dim>(fe_u, dim, fe_p, 1, FE_Nothing<dim>(fe_d.reference_cell()), dim));
    
    // Solid system: u (dim, deactivated), p (1, deactivated), d (dim)
    fe_collection.push_back(
      FESystem<dim>(FE_Nothing<dim>(fe_u.reference_cell()), dim, 
                   FE_Nothing<dim>(fe_p.reference_cell()), 1, fe_d, dim));

    quad_collection.push_back(QGaussSimplex<dim>(deg_u + 1)); // Stokes
    quad_collection.push_back(QGaussSimplex<dim>(deg_d + 1)); // Elasticity

    face_quad_collection.push_back(QGaussSimplex<dim - 1>(deg_u + 1));
    face_quad_collection.push_back(QGaussSimplex<dim - 1>(deg_d + 1));

    mapping_collection.push_back(MappingFE<dim>(FE_SimplexP<dim>(1)));
    mapping_collection.push_back(MappingFE<dim>(FE_SimplexP<dim>(1)));
  }

  // 2. DoF distribution
  dof_handler.reinit(triangulation);
  for (const auto &cell : dof_handler.active_cell_iterators())
  {
    if (cell->is_locally_owned())
    {
      if (cell->material_id() == DomainLabels::fluid)
        cell->set_active_fe_index(0);
      else
        cell->set_active_fe_index(1);
    }
  }
  dof_handler.distribute_dofs(fe_collection);

  locally_owned_dofs = dof_handler.locally_owned_dofs();
  locally_relevant_dofs = DoFTools::extract_locally_relevant_dofs(dof_handler);

  // 3. Constraints initialization
  {
    system_constraints.clear();
    system_constraints.reinit(locally_relevant_dofs);
    DoFTools::make_hanging_node_constraints(dof_handler, system_constraints);

    // Kinematic Coupling: u = 0 on interface
    std::vector<types::global_dof_index> face_dofs;
    for (const auto &cell : dof_handler.active_cell_iterators())
    {
      if (!cell->is_locally_owned() || cell->active_fe_index() != 0) continue;

      for (unsigned int f = 0; f < cell->n_faces(); ++f)
      {
        if (cell->face(f)->at_boundary()) continue;

        bool is_interface = false;
        if (cell->neighbor(f)->active_fe_index() == 1)
          is_interface = true;
        
        if (is_interface)
        {
          face_dofs.resize(cell->get_fe().n_dofs_per_face());
          cell->face(f)->get_dof_indices(face_dofs, 0);
          for (unsigned int i = 0; i < face_dofs.size(); ++i)
          {
            if (cell->get_fe().face_system_to_component_index(i).first < dim)
              system_constraints.add_line(face_dofs[i]);
          }
        }
      }
    }

    // Dirichlet Boundaries
    // Inlet
    InletFunction inlet_f;
    ComponentMask u_mask(dim + 1 + dim, false);
    for (unsigned int i = 0; i < dim; ++i) u_mask.set(i, true);
    VectorTools::interpolate_boundary_values(dof_handler, DomainLabels::inflow, inlet_f, system_constraints, u_mask);

    // Clamped Solid
    ComponentMask d_mask(dim + 1 + dim, false);
    for (unsigned int i = 0; i < dim; ++i) d_mask.set(dim + 1 + i, true);
    Functions::ZeroFunction<dim> zero_f(dim + 1 + dim);
    VectorTools::interpolate_boundary_values(dof_handler, DomainLabels::bottom_clamped, zero_f, system_constraints, d_mask);

    system_constraints.close();
  }

  // 4. Matrix and Vector allocation
  {
    DynamicSparsityPattern dsp(locally_relevant_dofs);
    
    Table<2, DoFTools::Coupling> coupling(dim + 1 + dim, dim + 1 + dim);
    for (unsigned int i = 0; i < dim + 1 + dim; ++i)
      for (unsigned int j = 0; j < dim + 1 + dim; ++j)
        coupling[i][j] = DoFTools::always;
    
    DoFTools::make_flux_sparsity_pattern(dof_handler, dsp, coupling, coupling);
    system_constraints.condense(dsp);
    SparsityTools::distribute_sparsity_pattern(dsp, locally_owned_dofs, MPI_COMM_WORLD, locally_relevant_dofs);

    system_matrix.reinit(locally_owned_dofs, locally_owned_dofs, dsp, MPI_COMM_WORLD);
    system_rhs.reinit(locally_owned_dofs, MPI_COMM_WORLD);
    solution_owned.reinit(locally_owned_dofs, MPI_COMM_WORLD);
    solution_relevant.reinit(locally_owned_dofs, locally_relevant_dofs, MPI_COMM_WORLD);
  }
}

void FSIProblem::assemble()
{
  pcout << "Assembling Monolithic System" << std::endl;
  system_matrix = 0.0;
  system_rhs = 0.0;

  hp::FEValues<dim> hp_fe_values(mapping_collection, fe_collection, quad_collection,
                                 update_values | update_gradients | update_quadrature_points | update_JxW_values);
  hp::FEFaceValues<dim> hp_fe_face_values(mapping_collection, fe_collection, face_quad_collection,
                                          update_values | update_gradients | update_normal_vectors | update_JxW_values);

  FEValuesExtractors::Vector velocity(0);
  FEValuesExtractors::Scalar pressure(dim);
  FEValuesExtractors::Vector displacement(dim + 1);

  std::vector<types::global_dof_index> dof_indices;

  for (const auto &cell : dof_handler.active_cell_iterators())
  {
    if (!cell->is_locally_owned()) continue;

    const unsigned int fe_index = cell->active_fe_index();
    const unsigned int dofs_per_cell = cell->get_fe().n_dofs_per_cell();
    FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell);
    Vector<double> cell_rhs(dofs_per_cell);

    hp_fe_values.reinit(cell);
    const FEValues<dim> &fe_values = hp_fe_values.get_present_fe_values();

    if (fe_index == 0) // Fluid
    {
      for (unsigned int q = 0; q < fe_values.n_quadrature_points; ++q)
      {
        for (unsigned int i = 0; i < dofs_per_cell; ++i)
        {
          for (unsigned int j = 0; j < dofs_per_cell; ++j)
          {
            cell_matrix(i, j) += (physical_params.viscosity * 
                                  scalar_product(fe_values[velocity].gradient(i, q), fe_values[velocity].gradient(j, q)) -
                                  fe_values[velocity].divergence(i, q) * fe_values[pressure].value(j, q) -
                                  fe_values[velocity].divergence(j, q) * fe_values[pressure].value(i, q)) *
                                 fe_values.JxW(q);
          }
        }
      }
    }
    else // Solid
    {
      for (unsigned int q = 0; q < fe_values.n_quadrature_points; ++q)
      {
        for (unsigned int i = 0; i < dofs_per_cell; ++i)
        {
          for (unsigned int j = 0; j < dofs_per_cell; ++j)
          {
            cell_matrix(i, j) += (physical_params.mu * 
                                  scalar_product(fe_values[displacement].gradient(i, q), fe_values[displacement].gradient(j, q)) +
                                  physical_params.lambda * 
                                  fe_values[displacement].divergence(i, q) * fe_values[displacement].divergence(j, q)) *
                                 fe_values.JxW(q);
          }
        }
      }
    }

    dof_indices.resize(dofs_per_cell);
    cell->get_dof_indices(dof_indices);
    system_constraints.distribute_local_to_global(cell_matrix, cell_rhs, dof_indices, system_matrix, system_rhs);

    // Interface Dynamic Coupling
    if (fe_index == 1) // On solid side
    {
      for (unsigned int f = 0; f < cell->n_faces(); ++f)
      {
        if (cell->face(f)->at_boundary()) continue;
        const auto neighbor = cell->neighbor(f);
        if (neighbor->active_fe_index() != 0) continue;

        hp_fe_face_values.reinit(cell, f); // Solid side
        const FEFaceValues<dim> &fe_face_solid = hp_fe_face_values.get_present_fe_values();
        
        hp_fe_face_values.reinit(neighbor, cell->neighbor_of_neighbor(f)); // Fluid side
        const FEFaceValues<dim> &fe_face_fluid = hp_fe_face_values.get_present_fe_values();

        FullMatrix<double> interface_matrix(dofs_per_cell, neighbor->get_fe().n_dofs_per_cell());
        
        for (unsigned int q = 0; q < fe_face_solid.n_quadrature_points; ++q)
        {
          const auto n = fe_face_fluid.normal_vector(q);
          for (unsigned int i = 0; i < dofs_per_cell; ++i)
          {
            const auto phi_d = fe_face_solid[displacement].value(i, q);
            for (unsigned int j = 0; j < neighbor->get_fe().n_dofs_per_cell(); ++j)
            {
              const auto traction = physical_params.viscosity * (fe_face_fluid[velocity].gradient(j, q) * n) -
                                    fe_face_fluid[pressure].value(j, q) * n;
              interface_matrix(i, j) += (traction * phi_d) * fe_face_solid.JxW(q);
            }
          }
        }
        
        std::vector<types::global_dof_index> neighbor_dof_indices(neighbor->get_fe().n_dofs_per_cell());
        neighbor->get_dof_indices(neighbor_dof_indices);
        system_constraints.distribute_local_to_global(interface_matrix, dof_indices, neighbor_dof_indices, system_matrix);
      }
    }
  }

  system_matrix.compress(VectorOperation::add);
  system_rhs.compress(VectorOperation::add);
}

void FSIProblem::solve()
{
  pcout << "Solving Linear System" << std::endl;
  SolverControl solver_control;
  TrilinosWrappers::SolverDirect direct_solver(solver_control);
  direct_solver.solve(system_matrix, solution_owned, system_rhs);
  system_constraints.distribute(solution_owned);
  solution_relevant = solution_owned;
}

void FSIProblem::output()
{
  pcout << "Writing Output Files" << std::endl;
  DataOut<dim> data_out;
  data_out.attach_dof_handler(dof_handler);

  std::vector<std::string> names;
  std::vector<DataComponentInterpretation::DataComponentInterpretation> interpretation;

  for (unsigned int d = 0; d < dim; ++d) { names.push_back("velocity"); interpretation.push_back(DataComponentInterpretation::component_is_part_of_vector); }
  names.push_back("pressure"); interpretation.push_back(DataComponentInterpretation::component_is_scalar);
  for (unsigned int d = 0; d < dim; ++d) { names.push_back("displacement"); interpretation.push_back(DataComponentInterpretation::component_is_part_of_vector); }

  data_out.add_data_vector(dof_handler, solution_relevant, names, interpretation);
  data_out.build_patches();

  const std::string filename = "output-fsi_" + Utilities::int_to_string(mpi_rank, 4) + ".vtu";
  std::ofstream out(filename);
  data_out.write_vtu(out);

  if (mpi_rank == 0)
  {
    std::vector<std::string> filenames;
    for (unsigned int i = 0; i < mpi_size; ++i) filenames.push_back("output-fsi_" + Utilities::int_to_string(i, 4) + ".vtu");
    std::ofstream master_out("output-fsi.pvtu");
    data_out.write_pvtu_record(master_out, filenames);
  }
}
