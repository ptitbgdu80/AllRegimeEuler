#include "AREuler.h"

int main()
{
  int nbr_elements_1D = 200;
  double t_final = 0.4;
  int choice_theta = 0;
  std::string file_name = "Riemann2Dtheta0";
  int choice_test_case = 1; // 0 = Vortex in a box, 1 = 2D-Riemann Problem, 2 = cas test Fanny B nul, 3 = cas 1D
  int choice_solver = 0; // 0 = prédiction correction, 1 = Rusanov


  Problem2D test(nbr_elements_1D, t_final, choice_theta, file_name, choice_test_case, choice_solver);

  test.Solve();

  return 0;
}
