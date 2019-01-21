#include "AREuler.h"

int main()
{
  int nbr_elements = 1000;
  double t_final = 3.1e-4;
  int choix_theta = 0;
  std::string file_name = "debug1D";

  Probleme1D test(nbr_elements, t_final, choix_theta, file_name);

  test.Solve();



  // int nbr_elements_1D = 50;
  // double t_final = 0.125;
  // int choix_theta = 2;
  // std::string file_name = "debug2D";
  //
  // Probleme2D test(nbr_elements_1D, t_final, choix_theta, file_name);
  //
  // test.Solve();

  return 0;
}
