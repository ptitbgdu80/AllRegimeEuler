#include "AREuler.h"

int main()
{
  int flow_choice = left_elem_flow;

  TroisVariables test(500, 0.0001, 2., flow_choice);

  test.ClassicGodunovMainLoop();

  return 0;
}
