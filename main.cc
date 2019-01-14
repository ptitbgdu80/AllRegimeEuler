#include "AREuler.h"

int main()
{
  int nbr_elements = 500;
  double delta_t = 0.001;
  double t_final = 2.;
  int flow_choice = left_elem_flow;
  std::string file_name = "resultats";

  TroisVariables test(nbr_elements, delta_t, t_final, flow_choice, file_name);

  test.ClassicFVMainLoop();

  return 0;
}
