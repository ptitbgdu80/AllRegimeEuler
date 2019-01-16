#include <vector>
#include <math.h>
#include <iostream>
#include <ostream>
#include <fstream>
#include <string>

class Probleme1D
{

protected:
  std::vector<std::vector<double>> _U, _Phi_interface;
  std::vector<double> _u, _u_star, _L, _Pi, _Pi_star, _left_bound_U, _right_bound_U, _a;
  double _nbr_elements, _delta_x, _Dt_on_Dx, _cfl, _time, _t_final, _left_bound_u, _right_bound_u, _left_bound_Pi, _right_bound_Pi;
  std::string _file_name;

public:
  Probleme1D(int nbr_elements, double t_final, std::string file_name);
  ~Probleme1D();
  void Update_u();
  void Update_Pi();
  void Update_a();
  void Update_u_star();
  void Update_Pi_star();
  void Update_L();
  void Update_Phi_interface();
  std::vector<double> LeftBoundValue();
  std::vector<double> RightBoundValue();
  void AcousticStep();
  void TransportStep();
  void SaveIteration(int time_it);
  void TimeIteration(int time_it);
  void Solve();
};

double PressureToRhoE(double p, double rho, double rho_u, double rho_v);
