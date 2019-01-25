#include <vector>
#include <math.h>
#include <iostream>
#include <ostream>
#include <fstream>
#include <string>

class Problem1D
{

protected:
  std::vector<std::vector<double>> _U, _Phi_interface;
  std::vector<double> _u, _u_star, _L, _Pi, _Pi_star, _left_bound_U, _right_bound_U, _a, _theta;
  double _nbr_elements, _delta_x, _Dt_on_Dx, _cfl, _time, _t_final, _left_bound_u, _right_bound_u, _left_bound_Pi, _right_bound_Pi, _down_bound_Pi, _up_bound_Pi;
  std::string _file_name;
  int _choice_theta;

public:
  Problem1D(int nbr_elements, double t_final, int choice_theta, std::string file_name);
  ~Problem1D();
  void Update_u();
  void Update_Pi();
  void Update_a();
  void Update_theta();
  void Update_u_star();
  void Update_Pi_star();
  void Update_L();
  void Update_Phi_interface();
  void Update_CL();
  void AcousticStep();
  void TransportStep();
  void SaveIteration(int time_it);
  void TimeIteration(int time_it);
  void Solve();
};

class Problem2D
{

protected:
  std::vector<std::vector<std::vector<double>>> _U, _Phi_interface_LR, _Phi_interface_DU, _flow_LR, _flow_DU, _Fx, _Fy;
  std::vector<std::vector<double>> _u, _v, _u_star, _v_star, _L, _Pi, _Pi_star_LR, _Pi_star_DU, _a_LR, _a_DU, _theta_LR, _theta_DU;
  std::vector<std::vector<double>> _left_bound_U, _right_bound_U, _down_bound_U, _up_bound_U, _max_c_LR, _max_c_DU, _left_bound_Fx, _right_bound_Fx, _down_bound_Fy, _up_bound_Fy;
  std::vector<double> _left_bound_u, _right_bound_u, _down_bound_v, _up_bound_v, _left_bound_Pi, _right_bound_Pi, _down_bound_Pi, _up_bound_Pi;
  double _nbr_elements_1D, _delta_s, _Dt_on_Ds, _cfl, _time, _t_final, _length;
  std::string _file_name;
  int _choice_theta, _choice_test_case, _choice_solver;

public:
  Problem2D(int nbr_elements_1D, double t_final, int choice_theta, std::string file_name, int choice_test_case, int choice_solver);
  ~Problem2D();
  void Update_u_v_Pi();
  void Update_a_u_v_star();
  void Update_theta();
  void Update_Pi_star();
  void Update_L();
  void Update_Phi_interface();
  void Update_CL();
  void AcousticStep();
  void TransportStep();
  void Update_F();
  void Update_max_c();
  void Update_flow();
  void RusanovStep();
  void SaveIteration(int time_it);
  void SaveCutsForRiemann2D();
  void TimeIteration(int time_it);
  void Solve();
};

double PressureToRhoE(double p, double rho, double rho_u, double rho_v);

double PressureToRhoE2D(double p, double rho, double rho_u, double rho_v);
