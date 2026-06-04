#include "FSIProblem.hpp"

// Main function.
int main(int argc, char *argv[])
{
  Utilities::MPI::MPI_InitFinalize mpi_init(argc, argv);

  const MPI_Comm comm = MPI_COMM_WORLD;
  const unsigned int rank = 
                Utilities::MPI::this_mpi_process(comm);

  const std::string mesh_file_name = "../mesh/mesh_fsi.msh";
  const unsigned int degree_velocity = 2;
  const unsigned int degree_pressure = 1;
  const unsigned int degree_displacement = 2;

  FSIProblem problem(mesh_file_name, degree_velocity, degree_pressure, degree_displacement);

  // Setup
  MPI_Barrier(comm);
  double t0 = MPI_Wtime();

  problem.setup();

  MPI_Barrier(comm);
  double t1 = MPI_Wtime();

  double setup_local = t1 - t0;
  double setup_time;

  MPI_Reduce(&setup_local, &setup_time, 1,
             MPI_DOUBLE, MPI_MAX, 0, comm);

  // Assemble
  MPI_Barrier(comm);
  t0 = MPI_Wtime();

  problem.assemble();

  MPI_Barrier(comm);
  t1 = MPI_Wtime();

  double assemble_local = t1 - t0;
  double assemble_time;

  MPI_Reduce(&assemble_local, &assemble_time, 1,
             MPI_DOUBLE, MPI_MAX, 0, comm);

  // Solve
  MPI_Barrier(comm);
  t0 = MPI_Wtime();

  problem.solve();

  MPI_Barrier(comm);
  t1 = MPI_Wtime();

  double solve_local = t1 - t0;
  double solve_time;

  MPI_Reduce(&solve_local, &solve_time, 1,
             MPI_DOUBLE, MPI_MAX, 0, comm);

  // Output
  problem.output();

  // Summary
  if (rank == 0)
  {
    std::cout << "Timing summary (max over MPI ranks)\n";
    std::cout << "setup    = " << setup_time    << " s\n";
    std::cout << "assemble = " << assemble_time << " s\n";
    std::cout << "solve    = " << solve_time    << " s\n";
    std::cout << "Total Time   = " << setup_time + 
                                      assemble_time + 
                                      solve_time  << " s\n";
    std::cout << "===============================================" << std::endl;
  }

  return 0;
}
