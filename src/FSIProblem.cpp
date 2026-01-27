#include "FSIProblem.hpp"


void
FSIProblem::setup()
{
  // Create the mesh.
  {
        std::cout << "Initializing the mesh" << std::endl;

    // Read the mesh from file.
    GridIn<dim> grid_in;
    grid_in.attach_triangulation(mesh);

    std::ifstream mesh_file(mesh_file_name);
    grid_in.read_msh(mesh_file);

    std::cout << "  Number of elements = " << mesh.n_active_cells()
              << std::endl;
  }

  pcout << "-----------------------------------------------" << std::endl;

  // Initialize the finite element space.
  {
    pcout << "Initializing the finite element space" << std::endl;

    // ----------------------------
    // Fluid: Stokes finite element
    // ----------------------------
    const FE_SimplexP<dim> fe_velocity(degree_velocity);
    const FE_SimplexP<dim> fe_pressure(degree_pressure);

    // Reference cell required to construct FE_Nothing
    const auto rc = fe_velocity.reference_cell();

    // FE_Nothing is used to deactivate unused components
    const FE_Nothing<dim> fe_nothing(rc); 

    // Stokes FE system:
    // [ u (dim components), p (1), d (dim components inactive) ]
    fe_stokes = std::make_unique<FESystem<dim>>(fe_velocity, dim, 
                                                fe_pressure, 1,
                                                fe_nothing, dim);

    // ------------------------------------
    // Solid: elasticity finite element
    // ------------------------------------
    const FE_SimplexP<dim> fe_displacement(degree_displacement);

    
    // Elasticity FE system:
    // [ u (dim components inactive), p (1 inactive), d (dim components) ]
    fe_elasticity = std::make_unique<FESystem<dim>>(fe_nothing, dim,
                                                    FE_Nothing<dim>(rc), 1,
                                                    fe_displacement, dim);
    pcout << "FE created" << std::endl;

    // ---------------------------------
    // hp::FECollection
    // index 0 -> fluid (Stokes)
    // index 1 -> solid (Elasticity)
    // ---------------------------------
    fe_collection.push_back(*fe_stokes);
    fe_collection.push_back(*fe_elasticity);

    pcout << "fe_collection size = " << fe_collection.size() << std::endl;

    // ===========================================================================
    // Mapping and quadrature
    // ===========================================================================
    mapping_collection.push_back(MappingFE<dim>(FE_SimplexP<dim>(1)));
    mapping_collection.push_back(MappingFE<dim>(FE_SimplexP<dim>(1)));

    quadrature_collection.push_back(QGaussSimplex<dim>(degree_velocity + 1));      // for stokes
    quadrature_collection.push_back(QGaussSimplex<dim>(degree_displacement + 1));  // for elasticity

    face_quadrature_collection.push_back(QGaussSimplex<dim-1>(degree_velocity + 1));
    face_quadrature_collection.push_back(QGaussSimplex<dim-1>(degree_displacement + 1));


    pcout << "  Velocity degree:           = " << fe_velocity.degree
          << std::endl;
    pcout << "  Pressure degree:           = " << fe_pressure.degree
          << std::endl;
    pcout << "  Displacement degree:       = " << fe_displacement.degree
          << std::endl;
   

    pcout << "  FE stokes components = " << fe_stokes->n_components() << std::endl;
    pcout << "  FE elasticity components = " << fe_elasticity->n_components() << std::endl;
  }

  pcout << "-----------------------------------------------" << std::endl;

  // Initialize the DoF handler.
  {
    pcout << "Initializing the DoF handler" << std::endl;

    pcout << "DoF: reinit..." << std::endl;

    //Attach DOFHandler to the mesh
    dof_handler.reinit(mesh);
    pcout << "DoF: set active_fe_index..." << std::endl;

    // Assign the active finite element based on material id
    // Fluid  -> FE index 0 (Stokes)
    // Solid  -> FE index 1 (Elasticity)
    for(const auto &cell : dof_handler.active_cell_iterators())
    {
      if(cell->material_id() == MaterialIds::Fluid)
        cell->set_active_fe_index(0); // stokes
      else if (cell->material_id() == MaterialIds::Solid)
        cell->set_active_fe_index(1); // elasticity
    }
pcout << "DoF: distribute..." << std::endl;

    // Distribute DoFs for the hp finite element space
    dof_handler.distribute_dofs(fe_collection);
pcout << "DoF: done distribute" << std::endl;
    /*
    // We want to reorder DoFs so that all velocity DoFs come first, and then
    // all pressure DoFs.
    std::vector<unsigned int> block_component(2*dim + 1, 0);
    for (unsigned int c=0; c<dim; ++c) block_component[c] = 0;     
    block_component[dim] = 1;                                       
    for (unsigned int c=dim+1; c<2*dim+1; ++c) block_component[c] = 2;
    DoFRenumbering::component_wise(dof_handler, block_component);
    */

    // In serial runs, locally owned and relevant DoFs coincide
    locally_owned_dofs = dof_handler.locally_owned_dofs();
    locally_relevant_dofs = DoFTools::extract_locally_relevant_dofs(dof_handler);

    pcout << "Total DoFs = " << dof_handler.n_dofs() << std::endl;

    // Besides the locally owned and locally relevant indices for the whole
    // system (velocity and pressure), we will also need those for the
    // individual velocity and pressure blocks.
    // std::vector<types::global_dof_index> dofs_per_block =
    //   DoFTools::count_dofs_per_fe_block(dof_handler, block_component);
/*
    block_owned_dofs.resize(3);
    block_relevant_dofs.resize(3);

    for (auto &is : block_owned_dofs)   is.clear();
    for (auto &is : block_relevant_dofs) is.clear();

    // Stokes
    ComponentMask mask_u_stokes(2*dim+1, false);
    for (unsigned int c=0;c<dim;++c) 
      mask_u_stokes.set(c,true);

    ComponentMask mask_p_stokes(2*dim+1, false);
    mask_p_stokes.set(dim,true);

    // Elasticity
    ComponentMask mask_d_solid(2*dim+1, false);
    for (unsigned int c=dim+1;c<2*dim+1;++c) 
      mask_d_solid.set(c,true);

    IndexSet u_owned(locally_owned_dofs.size()), p_owned(locally_owned_dofs.size()), d_owned(locally_owned_dofs.size());
    IndexSet u_rel(locally_relevant_dofs.size()), p_rel(locally_relevant_dofs.size()), d_rel(locally_relevant_dofs.size());

    u_owned = IndexSet(dof_handler.n_dofs());
    p_owned = IndexSet(dof_handler.n_dofs());
    d_owned = IndexSet(dof_handler.n_dofs());
    u_rel   = IndexSet(dof_handler.n_dofs());
    p_rel   = IndexSet(dof_handler.n_dofs());
    d_rel   = IndexSet(dof_handler.n_dofs());

    std::vector<types::global_dof_index> local_dof_indices;

    for (const auto &cell : dof_handler.active_cell_iterators())
    {
      local_dof_indices.resize(cell->get_fe().dofs_per_cell);
      cell->get_dof_indices(local_dof_indices);

      const unsigned int fe_index = cell->active_fe_index();

      if (fe_index == 0) // stokes cell
      {
        for (unsigned int i = 0; i < local_dof_indices.size(); ++i)
        {
          const auto gi = local_dof_indices[i];

          if (!locally_relevant_dofs.is_element(gi))
            continue;

          const unsigned int comp = cell->get_fe().system_to_component_index(i).first;

          if (comp < dim) u_rel.add_index(gi);
          else if (comp == dim) p_rel.add_index(gi);

          if (locally_owned_dofs.is_element(gi))
          {
            if (comp < dim) u_owned.add_index(gi);
            else if (comp == dim) p_owned.add_index(gi);
          }
        }
      }
      else if (fe_index == 1) // solid cell
      {
        for (unsigned int i = 0; i < local_dof_indices.size(); ++i)
        {
          const auto gi = local_dof_indices[i];

          if (!locally_relevant_dofs.is_element(gi))
            continue;

          const unsigned int comp = cell->get_fe().system_to_component_index(i).first;

          // here comp in [0..dim-1] corresponds to displacement components
          if (comp > dim && comp < 2*dim + 1) d_rel.add_index(gi);

          if (locally_owned_dofs.is_element(gi))
          {
            if (comp > dim && comp < 2*dim + 1) d_owned.add_index(gi);
          }
        }
      }
    }

    u_owned.compress(); p_owned.compress(); d_owned.compress();
    u_rel.compress();   p_rel.compress();   d_rel.compress();

    block_owned_dofs[0] = u_owned;
    block_owned_dofs[1] = p_owned;
    block_owned_dofs[2] = d_owned;

    block_relevant_dofs[0] = u_rel;
    block_relevant_dofs[1] = p_rel;
    block_relevant_dofs[2] = d_rel;




    pcout << "  Block DoFs: velocity = " << block_owned_dofs[0].n_elements()
          << " pressure = " << block_owned_dofs[1].n_elements()
          << " displacement = " << block_owned_dofs[2].n_elements() << std::endl;
*/
  }

  pcout << "-----------------------------------------------" << std::endl;

  // Initialize the linear system.
  {
    pcout << "Initializing the linear system" << std::endl;

    pcout << "  Initializing the sparsity pattern" << std::endl;


    // Velocity DoFs interact with other velocity DoFs (the weak formulation has
    // terms involving u times v), and pressure DoFs interact with velocity DoFs
    // (there are terms involving p times v or u times q). However, pressure
    // DoFs do not interact with other pressure DoFs (there are no terms
    // involving p times q). In Elasticity displacement DoFs interact with other displacement DoFs, 
    // and finally there is the interaction on the interface.
    // We build a table to store this information, so that
    // the sparsity pattern can be built accordingly.
    Table<2, DoFTools::Coupling> cell_coupling(2*dim + 1, 2*dim + 1);
    Table<2, DoFTools::Coupling> face_coupling(2*dim + 1, 2*dim + 1);

    for (unsigned int c = 0; c < 2*dim + 1; ++c)
      for (unsigned int d = 0; d < 2*dim + 1; ++d){
        cell_coupling[c][d] = DoFTools::none;
        face_coupling[c][d] = DoFTools::none;
      }
    // u-v
    for (unsigned int c = 0; c < dim; ++c)
      for (unsigned int d = 0; d < dim; ++d)
        cell_coupling[c][d] = DoFTools::always;

    // u-q e p-v
    for (unsigned int c = 0; c < dim; ++c)
    {
      cell_coupling[c][dim] = DoFTools::always;
      cell_coupling[dim][c] = DoFTools::always;
    }

    // d-b
    for (unsigned int c = dim+1; c < 2*dim+1; ++c)
      for (unsigned int d = dim+1; d < 2*dim+1; ++d)
        cell_coupling[c][d] = DoFTools::always;

    // Interface coupling
    for (unsigned int c = dim+1; c < 2*dim+1; ++c)
      for(unsigned int d = 0; d <= dim; ++d)
    {
      face_coupling[c][d] = DoFTools::always;
      //face_coupling[d][c] = DoFTools::always;
    }

    // velocity = 0 on the interface with constraints
    constraints.clear();
  
    //constraints.reinit(dof_handler.n_dofs());
  
    DoFTools::make_hanging_node_constraints(dof_handler, constraints);
  
    std::vector<types::global_dof_index> local_face_dof_indices(fe_stokes->n_dofs_per_face());
   

    for (const auto &cell : dof_handler.active_cell_iterators())
    {
      if (!cell->is_locally_owned())
        continue;

      if (cell->active_fe_index() != 0)
        continue;

      for (const auto face_no : cell->face_indices())
      {
        if (cell->face(face_no)->at_boundary())
          continue;

        bool face_is_on_interface = false;

        const auto neighbor = cell->neighbor(face_no);

        if (!neighbor->has_children() &&
            neighbor->active_fe_index() == 1)
        {
          face_is_on_interface = true;
        }
        else if (neighbor->has_children())
        {
          for (unsigned int sf = 0;
              sf < cell->face(face_no)->n_children();
              ++sf)
          {
            if (cell->neighbor_child_on_subface(face_no, sf)
                  ->active_fe_index() == 1)
            {
              face_is_on_interface = true;
              break;
            }
          }
        }

        if (face_is_on_interface)
        {
          // FE INDEX 0 = stokes
          cell->face(face_no)->get_dof_indices(local_face_dof_indices, 0);

          for (unsigned int i = 0; i < local_face_dof_indices.size(); ++i)
          {
            if (fe_stokes->face_system_to_component_index(i).first < dim)
            {
              constraints.add_line(local_face_dof_indices[i]);
              constraints.set_inhomogeneity(local_face_dof_indices[i], 0.0);

            }
          }
        }
      }

    }
    constraints.close();

    // ===========================================================================
    // Sparsity pattern and system matrices
    // ===========================================================================
    DynamicSparsityPattern dsp(dof_handler.n_dofs(), dof_handler.n_dofs());                         

    DoFTools::make_flux_sparsity_pattern(dof_handler, dsp, cell_coupling, face_coupling);

    constraints.condense(dsp);
    sparsity.copy_from(dsp);

    //TrilinosWrappers::BlockSparsityPattern sparsity;
    //sparsity.reinit(block_owned_dofs,   
    //                block_owned_dofs,   
    //                block_relevant_dofs,
    //                MPI_COMM_WORLD);
                    
    //sparsity.compress();
  
    /*
    Table<2, DoFTools::Coupling> coupling_pressure(2*dim + 1, 2*dim + 1);
    
    // We also build a sparsity pattern for the pressure mass matrix.
    for (unsigned int c = 0; c < 2*dim + 1; ++c)
      {
        for (unsigned int d = 0; d < 2*dim + 1; ++d)
          {
            if (c == dim && d == dim) // pressure-pressure term
              coupling_pressure[c][d] = DoFTools::always;
            else // other combinations
              coupling_pressure[c][d] = DoFTools::none;
          }
      }
    
    sparsity_pressure_mass.reinit(block_owned_dofs,   
                                  block_owned_dofs,   
                                  block_relevant_dofs);

    DoFTools::make_sparsity_pattern(dof_handler,
                                    coupling_pressure,
                                    sparsity_pressure_mass);
    sparsity_pressure_mass.compress();
*/
    pcout << "  Initializing the matrices" << std::endl;
    system_matrix.reinit(sparsity);
    //pressure_mass.reinit(sparsity_pressure_mass);

    pcout << "  Initializing the system right-hand side" << std::endl;
    system_rhs.reinit(dof_handler.n_dofs());
    pcout << "  Initializing the solution vector" << std::endl;
    solution_owned.reinit(dof_handler.n_dofs());
    solution.reinit(dof_handler.n_dofs());
  }
}

void
FSIProblem::assemble()
{
  pcout << "===============================================" << std::endl;
  pcout << "Assembling the system" << std::endl;

  // hp::FEValues allows handling different finite elements on different cells
  // (fluid vs solid) within the same loop.
  hp::FEValues<dim>     hp_fe_values(mapping_collection,
                          fe_collection,
                          quadrature_collection,
                          update_values | update_gradients |
                            update_quadrature_points | update_JxW_values);
  hp::FEFaceValues<dim> hp_fe_face_values(mapping_collection,
                                   fe_collection,
                                   face_quadrature_collection,
                                   update_values | update_gradients | 
                                   update_normal_vectors | update_JxW_values);
  
  
  // Global DoF indices on each cell
  std::vector<types::global_dof_index> dof_indices;

  system_matrix = 0.0;
  system_rhs    = 0.0;
  //pressure_mass = 0.0;

  // Component ordering:
  // [0 .. dim-1]     -> fluid velocity
  // [dim]            -> fluid pressure
  // [dim+1 .. 2*dim] -> solid displacement
  FEValuesExtractors::Vector velocity(0);
  FEValuesExtractors::Scalar pressure(dim);
  FEValuesExtractors::Vector displacement(dim+1);

  for (const auto &cell : dof_handler.active_cell_iterators())
    {
      //if (!cell->is_locally_owned())
        //continue;

      const unsigned int fe_index = cell->active_fe_index();
      const unsigned int dofs_per_cell = cell->get_fe().n_dofs_per_cell();

      FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell);
      //FullMatrix<double> cell_pressure_mass_matrix(dofs_per_cell, dofs_per_cell);
      Vector<double>     cell_rhs(dofs_per_cell);

      cell_matrix               = 0.0;
      cell_rhs                  = 0.0;
      //cell_pressure_mass_matrix = 0.0;

      dof_indices.resize(dofs_per_cell);

      // Reinitialize FEValues for the current cell
      hp_fe_values.reinit(cell);
      const FEValues<dim> &fe_values = hp_fe_values.get_present_fe_values();

          // FLUID CELL: Stokes equations
          if(fe_index == 0)
          {
            for(unsigned int q = 0; q < fe_values.n_quadrature_points; ++q)
            {
              for (unsigned int i = 0; i < dofs_per_cell; ++i)
                {
                  for (unsigned int j = 0; j < dofs_per_cell; ++j)
                    {
                    // Viscosity term
                    cell_matrix(i, j) +=
                      nu *
                      scalar_product(fe_values[velocity].gradient(i, q),
                                    fe_values[velocity].gradient(j, q)) *
                      fe_values.JxW(q);

                      // Pressure term in the momentum equation
                    cell_matrix(i, j) -= fe_values[velocity].divergence(i, q) *
                                        fe_values[pressure].value(j, q) *
                                        fe_values.JxW(q);

                    // Pressure term in the continuity equation
                    cell_matrix(i, j) -= fe_values[velocity].divergence(j, q) *
                                        fe_values[pressure].value(i, q) *
                                        fe_values.JxW(q);

                    // Pressure mass matrix.
                    //cell_pressure_mass_matrix(i, j) +=
                    //(fe_values[pressure].value(i, q) *
                    //fe_values[pressure].value(j, q) / nu) * fe_values.JxW(q);
                    
                    }
                }
            }
            /* Neumann boundary condition on outflow
            for(unsigned int f = 0; f < cell->n_faces(); ++f)
            {
              if(!cell->face(f)->at_boundary())
                continue;

              if(cell->face(f)->boundary_id() != BoundaryIds::Outflow)
                continue;

              hp_fe_face_values.reinit(cell, f);
              const FEFaceValues<dim> &fe_face = hp_fe_face_values.get_present_fe_values();

              for(unsigned int q = 0; q < fe_face.n_quadrature_points; ++q)
              {
                for(unsigned int i = 0; i < dofs_per_cell; ++i)
                {
                  cell_rhs(i) +=
                                p_out * 
                                scalar_product(fe_face.normal_vector(q), 
                                fe_face[velocity].value(i,q)) *
                                fe_face.JxW(q);
                }
              }
            }*/
          }

          // SOLID CELL: Linear elasticity
          if(fe_index == 1)
          {
            for(unsigned int q = 0; q < fe_values.n_quadrature_points; ++q)
            {
              for (unsigned int i = 0; i < dofs_per_cell; ++i)
                {
                  for (unsigned int j = 0; j < dofs_per_cell; ++j)
                    {
                    cell_matrix(i, j) +=
                                        (mu *
                                        scalar_product(fe_values[displacement].gradient(j, q),
                                        fe_values[displacement].gradient(i, q)) +
                                        lambda * fe_values[displacement].divergence(j, q) *
                                        fe_values[displacement].divergence(i, q)) *
                                        fe_values.JxW(q);
                    }
                    // No forcing term
                }
            }
          }
      
      //dof_indices.resize(cell->get_fe().n_dofs_per_cell());
      cell->get_dof_indices(dof_indices);
      //system_matrix.add(dof_indices, cell_matrix);
      //system_rhs.add(dof_indices, cell_rhs);
      constraints.distribute_local_to_global(cell_matrix,
                                       cell_rhs,
                                       dof_indices,
                                       system_matrix,
                                       system_rhs);

      //if(fe_index == 0)
        //  pressure_mass.add(dof_indices, cell_pressure_mass_matrix);

      // FLUID–SOLID INTERFACE COUPLING
      if(fe_index == 1)
      {
        for(unsigned int f = 0; f < cell->n_faces(); ++f)
        {
          if(cell->face(f)->at_boundary())
            continue;
          
          const auto neighbor = cell->neighbor(f);
          if(neighbor->active_fe_index() != 0)
            continue;

          // FEValues on solid side  
          hp_fe_face_values.reinit(cell,f);
          const FEFaceValues<dim> &fe_face_solid = hp_fe_face_values.get_present_fe_values();

          // FEValues on fluid side
          const unsigned int f_fluid = cell->neighbor_of_neighbor(f);
          hp_fe_face_values.reinit(neighbor, f_fluid);
          const FEFaceValues<dim> &fe_face_fluid = hp_fe_face_values.get_present_fe_values();

          const unsigned int dofs_solid = cell->get_fe().dofs_per_cell;
          const unsigned int dofs_fluid = neighbor->get_fe().dofs_per_cell;

          FullMatrix<double> local_interface_matrix(dofs_solid, dofs_fluid);
          FullMatrix<double> local_interface_matrix_T(dofs_fluid, dofs_solid);
          local_interface_matrix = 0.0;
          local_interface_matrix_T = 0.0;

          for(unsigned int q = 0; q < fe_face_solid.n_quadrature_points; ++q)
          {
            const Tensor<1, dim> n_fluid = fe_face_fluid.normal_vector(q);

            for(unsigned int i = 0; i < dofs_solid; ++i)
            {
              const Tensor<1, dim> phi_d_i = fe_face_solid[displacement].value(i,q);

              for(unsigned int j = 0; j < dofs_fluid; ++j)
              {
                const Tensor<1, dim> t_f = 
                                          (nu * (fe_face_fluid[velocity].gradient(j,q) * n_fluid)
                                          - fe_face_fluid[pressure].value(j,q) * n_fluid);

                local_interface_matrix(i,j) +=
                                              (t_f * phi_d_i) * fe_face_solid.JxW(q);

                local_interface_matrix_T(j,i) +=
                                              (t_f * phi_d_i) * fe_face_solid.JxW(q);
              }
            }
          }

          std::vector<types::global_dof_index> dofs_solid_indices(dofs_solid), dofs_fluid_indices(dofs_fluid);
          cell->get_dof_indices(dofs_solid_indices);
          neighbor->get_dof_indices(dofs_fluid_indices);

          constraints.distribute_local_to_global(local_interface_matrix,
                                                dofs_solid_indices,
                                                dofs_fluid_indices,
                                                system_matrix);

          //constraints.distribute_local_to_global(local_interface_matrix_T,
            //                                    dofs_fluid_indices,
              //                                  dofs_solid_indices,
                //                                system_matrix);

        }
      }
    }
      

  //system_matrix.compress(VectorOperation::add);
  //system_rhs.compress(VectorOperation::add);
  //pressure_mass.compress(VectorOperation::add);


  // Dirichlet boundary conditions.
  {
    std::map<types::global_dof_index, double>           boundary_values;
    std::map<types::boundary_id, const Function<dim> *> boundary_functions;

    // Fluid velocity BCs
    ComponentMask mask_velocity_stokes(fe_stokes->n_components(), false);
    for(unsigned int c = 0; c < dim; ++c)
      mask_velocity_stokes.set(c, true);

    boundary_functions[BoundaryIds::Inflow] = &inlet_velocity;

    VectorTools::interpolate_boundary_values(dof_handler,
                                             boundary_functions,
                                             boundary_values,
                                             mask_velocity_stokes);

    //boundary_functions.clear();
    //Functions::ZeroFunction<dim> zero_function_stokes(fe_stokes->n_components());

    //boundary_functions[BoundaryIds::LeftWallFluid] = &zero_function_stokes;
    //boundary_functions[BoundaryIds::RightWallFluid] = &zero_function_stokes;
    //VectorTools::interpolate_boundary_values(dof_handler,
    //                                        boundary_functions,
    //                                         boundary_values,
    //                                         mask_velocity_stokes);

    MatrixTools::apply_boundary_values(
                                      boundary_values, 
                                      system_matrix, 
                                      solution_owned, 
                                      system_rhs, 
                                      false);

    boundary_functions.clear();

    // Solid displacement BCs                                  
    ComponentMask mask_displacement_elasticity(fe_elasticity->n_components(), false);
    for (unsigned int c=dim+1; c<2*dim+1; ++c) 
      mask_displacement_elasticity.set(c,true);

    Functions::ZeroFunction<dim> zero_function_elasticity(fe_elasticity->n_components());

    //boundary_functions[BoundaryIds::LeftWallSolid] = &zero_function_elasticity;
    //boundary_functions[BoundaryIds::RightWallSolid] = &zero_function_elasticity;
    boundary_functions[BoundaryIds::BottomWall] = &zero_function_elasticity;

    VectorTools::interpolate_boundary_values(dof_handler,
                                             boundary_functions,
                                             boundary_values,
                                             mask_displacement_elasticity);

    MatrixTools::apply_boundary_values(
                                      boundary_values,
                                      system_matrix, 
                                      solution_owned, 
                                      system_rhs, 
                                      false);
  }
}

void
FSIProblem::solve()
{
  pcout << "===============================================" << std::endl;

  //SolverControl solver_control(9000, 1e-6 * system_rhs.l2_norm());

  //SolverGMRES<Vector<double>> solver(solver_control);

  //PreconditionBlockDiagonal preconditioner;
  //preconditioner.initialize(system_matrix.block(0, 0),
  //                          pressure_mass.block(1, 1),
  //                          system_matrix.block(2,2));

  //PreconditionIdentity preconditioner;
  

  pcout << "Solving the linear system" << std::endl;

  // Serial implementation with direct solver
  SparseDirectUMFPACK direct_solver;
  direct_solver.initialize(system_matrix);
  direct_solver.vmult(solution_owned, system_rhs);
  //solver.solve(system_matrix, solution_owned, system_rhs, preconditioner);
  //pcout << "  " << solver_control.last_step() << " GMRES iterations"
       // << std::endl;

  // Distribute constrained degrees of freedom
  constraints.distribute(solution_owned);

  // Copy solution to the final solution vector
  solution = solution_owned;
  //solution.update_ghost_values();
}

void
FSIProblem::output()
{
  pcout << "===============================================" << std::endl;

  DataOut<dim> data_out;

  data_out.attach_dof_handler(dof_handler);

   std::vector<std::string> names;
  std::vector<DataComponentInterpretation::DataComponentInterpretation> interpretation;

  // velocity (dim components)
  for (unsigned int d = 0; d < dim; ++d)
  {
    names.push_back("velocity");
    interpretation.push_back(DataComponentInterpretation::component_is_part_of_vector);
  }

  // pressure (1 component)
  names.push_back("pressure");
  interpretation.push_back(DataComponentInterpretation::component_is_scalar);

  // displacement (dim components)
  for (unsigned int d = 0; d < dim; ++d)
  {
    names.push_back("displacement");
    interpretation.push_back(DataComponentInterpretation::component_is_part_of_vector);
  }

  data_out.add_data_vector(dof_handler,
                           solution,
                           names,
                           interpretation);
/*
  std::vector<DataComponentInterpretation::DataComponentInterpretation>
    interpretation(dim,
                   DataComponentInterpretation::component_is_part_of_vector);
  interpretation.push_back(DataComponentInterpretation::component_is_scalar);
  for (unsigned int d = 0; d < dim; ++d)
  interpretation.push_back(DataComponentInterpretation::component_is_part_of_vector);

  std::vector<std::string> names(dim, "velocity");
  names.push_back("pressure");
  for(unsigned d = 0; d < dim; ++d)
  { 
    names.push_back("displacement");
  }

  data_out.add_data_vector(dof_handler, solution, names, interpretation);
*/

  std::vector<unsigned int> partition_int(mesh.n_active_cells());
  GridTools::get_subdomain_association(mesh, partition_int);
  const Vector<double> partitioning(partition_int.begin(), partition_int.end());
  data_out.add_data_vector(partitioning, "partitioning");

  data_out.build_patches();

  const std::string output_file_name = "output-FSIProblem";
  std::ofstream out(output_file_name + ".vtu");
  data_out.write_vtu(out);

  pcout << "Output written to " << output_file_name << std::endl;
  pcout << "===============================================" << std::endl;
}