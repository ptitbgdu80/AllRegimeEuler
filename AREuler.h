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
  std::vector<double> _u, _u_star, _L, _Pi, _Pi_star, _left_bound_U, _right_bound_U, _a, _theta;
  double _nbr_elements, _delta_x, _Dt_on_Dx, _cfl, _time, _t_final, _left_bound_u, _right_bound_u, _left_bound_Pi, _right_bound_Pi, _down_bound_Pi, _up_bound_Pi;
  std::string _file_name;
  int _choix_theta;

public:
  Probleme1D(int nbr_elements, double t_final, int choix_theta, std::string file_name);
  ~Probleme1D();
  void Update_u();
  void Update_Pi();
  void Update_a();
  void Update_theta();
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

class Probleme2D
{

protected:
  std::vector<std::vector<std::vector<double>>> _U, _Phi_interface_LR, _Phi_interface_DU;
  std::vector<std::vector<double>> _u, _v, _u_star, _v_star, _L, _Pi, _Pi_star_LR, _Pi_star_DU, _a_LR, _a_DU, _theta_LR, _theta_DU, _left_bound_U, _right_bound_U, _down_bound_U, _up_bound_U;
  std::vector<double> _left_bound_u, _right_bound_u, _down_bound_u, _up_bound_u, _left_bound_v, right_bound_v, _down_bound_v, _up_bound_v, _left_bound_Pi, _right_bound_Pi, _down_bound_Pi, _up_bound_Pi;
  double _nbr_elements_1D, _delta_s, _Dt_on_Ds, _cfl, _time, _t_final;
  std::string _file_name;
  int _choix_theta;

public:
  Probleme2D(int nbr_elements_1D, double t_final, int choix_theta, std::string file_name);
  ~Probleme2D();
  void Update_u_v_Pi();
  void Update_a_u_v_star();
  void Update_theta();
  void Update_Pi_star();
  void Update_L();
  void Update_Phi_interface();
  // std::vector<double> LeftBoundValue();
  // std::vector<double> RightBoundValue();
  void AcousticStep();
  // void TransportStep();
  // void SaveIteration(int time_it);
  // void TimeIteration(int time_it);
  // void Solve();
};

double PressureToRhoE(double p, double rho, double rho_u, double rho_v);

double PressureToRhoE2D(double p, double rho, double rho_u, double rho_v);
