#include "AREuler.h"

int main()
{
  int nbr_elements = 1000;
  double t_final = 3.1e-4;
  std::string file_name = "resultats";

  Probleme1D test(nbr_elements, t_final, file_name);

  test.Solve();

  return 0;
}
