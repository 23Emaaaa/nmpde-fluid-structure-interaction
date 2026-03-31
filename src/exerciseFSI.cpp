#include <deal.II/base/mpi.h>  // needed to run in parallel with MPI

#include "FSIProblem.hpp"

// Main function.
int main(int argc, char* argv[]) {
    // This wakes up the MPI network, assigns the ranks, and safely tears
    // down the network when the program finishes.
    dealii::Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);

    const std::string mesh_file_name = "../mesh/mesh_fsi.msh";
    const unsigned int degree_velocity = 2;
    const unsigned int degree_pressure = 1;
    const unsigned int degree_displacement = 2;

    FSIProblem problem(mesh_file_name, degree_velocity, degree_pressure,
                       degree_displacement);

    problem.setup();
    problem.assemble();
    problem.solve();
    problem.output();

    return 0;
}