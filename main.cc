#include "AREuler.h"

int main()
{
  // int nbr_elements = 1000;
  // double t_final = 3.1e-4;
  // int choice_theta = 0;
  // std::string file_name = "debug1D";
  //
  // Problem1D test(nbr_elements, t_final, choice_theta, file_name);
  //
  // test.Solve();


  int nbr_elements_1D = 300;
  double t_final = 3.1e-4;
  int choice_theta = 0;
  std::string file_name = "test1D";
  int choice_test_case = 3; // 0 = Vortex in a box, 1 = 2D-Riemann Problem, 2 = cas test Fanny B nul, 3 = cas 1D
  int choice_solver = 0; // 0 = prédiction correction, 1 = Rusanov


  Problem2D test(nbr_elements_1D, t_final, choice_theta, file_name, choice_test_case, choice_solver);

  test.Solve();

  return 0;
}
