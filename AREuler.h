#include <vector>
#include <math.h>
#include <iostream>
#include <string>

class Probleme1D
{

protected:
  std::vector<std::vector<double>> _U, _W;
  double _nbr_elements, _delta_x, _delta_t, _t_final;
  std::string _file_name;

public:
  Probleme1D(int nbr_elements, double delta_t, double t_final, std::string file_name);
  ~Probleme1D();
  std::vector<double> RightInterfaceFlow(int elem_i);
  std::vector<double> LeftBoundFlow();
  std::vector<double> RightBoundFlow();
  void AcousticStep();
  void ClassicFVMainLoop();
  void SaveIteration(int time_it);

};
