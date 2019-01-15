#include <vector>
#include <math.h>
#include <iostream>
#include <string>

class Probleme1D
{

protected:
  std::vector<std::vector<double>> _U;
  std::vector<double> _u, _u_star, _L, _Pi, _Pi_star, _left_bound_U, _right_bound_U, _a;
  double _nbr_elements, _Dt_on_Dx, _t_final, _left_bound_u, _right_bound_u, _left_bound_Pi, _right_bound_Pi;
  std::string _file_name;

public:
  Probleme1D(int nbr_elements, double delta_t, double t_final, std::string file_name);
  ~Probleme1D();
  void Update_u();
  void Update_Pi();
  void Update_u_star();
  void Update_Pi_star();
  void Update_L();
  std::vector<double> LeftBoundValue();
  std::vector<double> RightBoundValue();
  void AcousticStep();
  void SaveIteration(int time_it);

};
