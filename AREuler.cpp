#include "AREuler.h"

//Constructeur
Probleme1D::Probleme1D(int nbr_elements, double delta_t, double t_final, std::string file_name)
{
  _nbr_elements = nbr_elements;
  _t_final = t_final;
  _file_name = file_name;
  _Dt_on_Dx = delta_t*nbr_elements;

  _U.resize(nbr_elements); // U contient les variables conservatives (rho, rho*u, rho*v, rho*E)
  _u.resize(nbr_elements);
  _a.resize(nbr_elements + 1);
  _u_star.resize(nbr_elements + 1); //_u_star[j] = u*_(j-1/2) (u* sur l'interface gauche de l'élément contenant xj)
  _L.resize(nbr_elements); // L[j] = 1 + Dt/Dx*(u*_(j+1/2) - u*_(j-1/2))
  _Pi.resize(nbr_elements); // Pi contient les pressions calculées à partir de U
  _Pi_star.resize(nbr_elements + 1); //_Pi_star[j] = Pi*_(j-1/2) (Pi* sur l'interface gauche de l'élément contenant xj)

  _left_bound_Pi = 1.0;
  _right_bound_Pi = 1.0;

  double x = test(1.0,1.0,1.0,1.0);
  x = x + 1;

  _left_bound_U = LeftBoundValue();
  _right_bound_U = RightBoundValue();

  _left_bound_u = _left_bound_U[1]/_left_bound_U[0];
  _right_bound_u = _right_bound_U[1]/_right_bound_U[0];

  for (int elem_j = 0; elem_j < nbr_elements; elem_j++)
  {
    _U[elem_j].resize(4);
    _U[elem_j][0] = 1.0;
    _U[elem_j][1] = 10.0;
    _U[elem_j][2] = 5.0;
    _U[elem_j][3] = PressureToRhoE(1.0, _U[elem_j][0], _U[elem_j][1], _U[elem_j][2]);
  }
}



//Destructeur
Probleme1D::~Probleme1D()
{}



//Methodes
double PressureToRhoE(double p, double rho, double rho_u, double rho_v)
{
  double rho_E = 2.5*p - 0.5*(rho_u*rho_u + rho_v*rho_v)/rho;
  return rho_E;
}

double test(double p, double rho, double rho_u, double rho_v)
{
  return 1.0;
}

void Probleme1D::Update_u()
{
  for (int elem_j = 0; elem_j < _nbr_elements; elem_j++)
  {
    _u[elem_j] = _U[elem_j][1]/_U[elem_j][0];
  }
}

void Probleme1D::Update_Pi()
{
  for (int elem_j = 0; elem_j < _nbr_elements; elem_j++)
  {
    _Pi[elem_j] = 0.4*(_U[elem_j][3] - (_U[elem_j][1]*_U[elem_j][1] + _U[elem_j][2]*_U[elem_j][2])/(2*_U[elem_j][0])); // p = (gamma - 1)*rho*e avec gamma = 1.4 et e = E - |u²|
    if (_Pi[elem_j] < 0)
    {
      std::cout << "Attention pression négative pour l'élément " << elem_j << std::endl;
    }
  }
}

void Probleme1D::Update_a()
{
  // a = 1.5*max(rho_L*c_L, rho_R*c_R) , c = sqrt(gamma*p/rho) => rho*c = sqrt(gamma*p*rho)
  double rhoLcL = sqrt(1.4*_left_bound_Pi*_left_bound_U[0]);
  double rhoRcR = sqrt(1.4*_Pi[0]*_U[0][0]);
  _a[0] = 1.5*std::max(rhoLcL, rhoRcR);

  std::cout << "P_gauche = " << _left_bound_Pi << std::endl;
  std::cout << "rho_gauche =  " << _left_bound_U[0] << std::endl;
  std::cout << "a[0] = " << _a[0] << std::endl;

  for (int elem_j = 1; elem_j < _nbr_elements; elem_j++)
  {
    rhoLcL = rhoRcR;
    rhoRcR = sqrt(1.4*_Pi[elem_j]*_U[elem_j][0]);
    _a[elem_j] = 1.5*std::max(rhoLcL, rhoRcR);
  }

  rhoLcL = rhoRcR;
  rhoRcR = sqrt(1.4*_right_bound_Pi*_right_bound_U[0]);
  _a[_nbr_elements] = 1.5*std::max(rhoLcL, rhoRcR);
}

void Probleme1D::Update_u_star()
{
  _u_star[0] = 0.5*(_u[0] + _left_bound_u - (_Pi[0] - _left_bound_Pi)/_a[0]);
  for (int elem_j = 1; elem_j < _nbr_elements; elem_j++)
  {
    _u_star[elem_j] = 0.5*(_u[elem_j] + _u[elem_j-1] - (_Pi[elem_j] - _Pi[elem_j-1])/_a[elem_j]);
  }
  _u_star[_nbr_elements] = 0.5*(_right_bound_u + _u[_nbr_elements-1] - (_right_bound_Pi - _Pi[_nbr_elements-1])/_a[_nbr_elements]);
}

void Probleme1D::Update_Pi_star()
{
  double a = 1.0;
  _Pi_star[0] = 0.5*(_Pi[0] + _left_bound_Pi - a*(_u[0] - _left_bound_u));
  for (int elem_j = 1; elem_j < _nbr_elements; elem_j++)
  {
    _Pi_star[elem_j] = 0.5*(_Pi[elem_j] + _Pi[elem_j-1] - a*(_u[elem_j] - _u[elem_j-1]));
  }
  _Pi_star[_nbr_elements] = 0.5*(_right_bound_Pi + _Pi[_nbr_elements-1] - a*(_right_bound_u - _u[_nbr_elements-1]));
}

void Probleme1D::Update_L()
{
  for (int elem_j = 0; elem_j < _nbr_elements; elem_j++)
  {
    _L[elem_j] = 1. + _Dt_on_Dx*(_u_star[elem_j+1] - _u_star[elem_j]);
  }
}

std::vector<double> Probleme1D::LeftBoundValue()
{
  double rho, u, v, rho_u, rho_v, rho_E;
  rho = 1.0;
  u = 10.0;
  v = 5.0;

  rho_u = rho*u;
  rho_v = rho*v;
  rho_E = PressureToRhoE(_left_bound_Pi, rho, rho_u, rho_v);
  return {rho, rho_u, rho_v, rho_E};
}

std::vector<double> Probleme1D::RightBoundValue()
{
  double rho, u, v, rho_u, rho_v, rho_E;
  rho = 1.0;
  u = 10.0;
  v = 5.0;

  rho_u = rho*u;
  rho_v = rho*v;
  rho_E = PressureToRhoE(_right_bound_Pi, rho, rho_u, rho_v);
  return {rho, rho_u, rho_v, rho_E};
}

void Probleme1D::AcousticStep()
{
  std::cout << "debut de l'étape acoustique, u = " << std::endl;
  for (int elem_j = 0; elem_j < _nbr_elements; elem_j++)
  {
    std::cout << _U[elem_j][1]/_U[elem_j][0] << std::endl;
  }

  Update_u();
  Update_Pi();
  Update_a();
  Update_u_star();
  Update_Pi_star();
  Update_L();

  for (int elem_j = 0; elem_j < _nbr_elements; elem_j++)
  {
    _U[elem_j][0] = _U[elem_j][0]/_L[elem_j];
    _U[elem_j][1] = (_U[elem_j][1] - _Dt_on_Dx*(_Pi_star[elem_j+1] - _Pi_star[elem_j]))/_L[elem_j];
    _U[elem_j][2] = _U[elem_j][2]/_L[elem_j];
    _U[elem_j][3] = (_U[elem_j][3] - _Dt_on_Dx*(_Pi_star[elem_j+1]*_u_star[elem_j+1] - _Pi_star[elem_j]*_u_star[elem_j]))/_L[elem_j];
  }

  std::cout << "fin de l'étape acoustique, u = " << std::endl;
  for (int elem_j = 0; elem_j < _nbr_elements; elem_j++)
  {
    std::cout << _U[elem_j][1]/_U[elem_j][0] << std::endl;
  }

}

void Probleme1D::SaveIteration(int time_it)
{}
