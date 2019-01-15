#include "AREuler.h"

int main()
{
  int nbr_elements = 10;
  double delta_t = 0.1;
  double t_final = 2.;
  std::string file_name = "resultats";

  Probleme1D test(nbr_elements, delta_t, t_final, file_name);

  test.AcousticStep();

  return 0;
}
