#include "FSIProblem.hpp"

// Main function.
int main()
{
  const std::string  mesh_file_name  = "../mesh/mesh_fsi.msh";
  const unsigned int degree_velocity = 2;
  const unsigned int degree_pressure = 1;
  const unsigned int degree_displacement = 2;

  FSIProblem problem(mesh_file_name, degree_velocity, degree_pressure, degree_displacement);

  problem.setup();
  problem.assemble();
  problem.solve();
  problem.output();

  return 0;
}