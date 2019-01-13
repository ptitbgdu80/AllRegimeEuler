#include "AREuler.h"

int main()
{
  int choix_flux = left_elem_flow;

  TroisVariables test(50, 0.01, choix_flux);

  test.ClassicGodunovIteration();

  return 0;
}
