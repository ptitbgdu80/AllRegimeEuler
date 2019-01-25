#include "AREuler.h"

int main()
{
  int nbr_elements_1D = 100;
  double t_final = 0.125;
  int choice_theta = 1;
  std::string file_name = "Vortextheta1PC";
  int choice_test_case = 0; // 0 = Vortex in a box, 1 = 2D-Riemann Problem, 2 = cas test Fanny B nul, 3 = cas 1D
  int choice_solver = 0; // 0 = prédiction correction, 1 = Rusanov


  Problem2D test(nbr_elements_1D, t_final, choice_theta, file_name, choice_test_case, choice_solver);

  test.Solve();

  return 0;
}
