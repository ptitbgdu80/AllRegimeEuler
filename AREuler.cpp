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
  _u_star.resize(nbr_elements + 1); //_u_star[j] = u*_(j-1/2) (u* sur l'interface gauche de l'élément contenant xj)
  _L.resize(nbr_elements); // L[j] = 1 + Dt/Dx*(u*_(j+1/2) - u*_(j-1/2))
  _Pi.resize(nbr_elements); // Pi contient les pressions calculées à partir de U
  _Pi_star.resize(nbr_elements + 1); //_Pi_star[j] = Pi*_(j-1/2) (Pi* sur l'interface gauche de l'élément contenant xj)

  _left_bound_U = LeftBoundValue();
  _right_bound_U = RightBoundValue();

  _left_bound_u = _left_bound_U[1]/_left_bound_U[0];
  _right_bound_u = _right_bound_U[1]/_right_bound_U[0];
  _left_bound_Pi = 0.4*(_left_bound_U[3] - (_left_bound_U[1]*_left_bound_U[1] + _left_bound_U[2]*_left_bound_U[2])/_left_bound_U[0]);
  _right_bound_Pi = 0.4*(_right_bound_U[3] - (_right_bound_U[1]*_right_bound_U[1] + _right_bound_U[2]*_right_bound_U[2])/_right_bound_U[0]);

  for (int elem_j = 0; elem_j < nbr_elements; elem_j++)
  {
    _U[elem_j].resize(4);
    for (int jVar = 0; jVar < 5; jVar++)
    {
      _U[elem_j][jVar] = 10;
    }
  }
}



//Destructeur
Probleme1D::~Probleme1D()
{}



//Methodes
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
    _Pi[elem_j] = 0.4*(_U[elem_j][3] - (_U[elem_j][1]*_U[elem_j][1] + _U[elem_j][2]*_U[elem_j][2])/_U[elem_j][0]); // p = (gamma - 1)*rho*e avec gamma = 1.4 et e = E - |u²|
  }
}

void Probleme1D::Update_u_star()
{
  double a = 1.;
  _u_star[0] = 0.5*(_u[0] + _left_bound_u - (_Pi[0] - _left_bound_Pi)/a);
  for (int elem_j = 1; elem_j < _nbr_elements; elem_j++)
  {
    _u_star[elem_j] = 0.5*(_u[elem_j] + _u[elem_j-1] - (_Pi[elem_j] - _Pi[elem_j-1])/a);
  }
  _u_star[_nbr_elements] = 0.5*(_right_bound_u + _u[_nbr_elements-1] - (_right_bound_Pi - _Pi[_nbr_elements-1])/a);
}

void Probleme1D::Update_Pi_star()
{
  double a = 1.;
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
  double rho, u, v, E;
  rho = 1.2; //masse volumique de l'air
  u = 10.;
  v = 5.;
  E = 1.;
  return {rho, rho*u, rho*v, rho*E};
}

std::vector<double> Probleme1D::RightBoundValue()
{
  double rho, u, v, E;
  rho = 1.2; //masse volumique de l'air
  u = 10.;
  v = 5.;
  E = 1.;
  return {rho, rho*u, rho*v, rho*E};
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
  Update_u_star();
  Update_Pi_star();
  Update_L();

  for (int elem_j = 0; elem_j < _nbr_elements; elem_j++)
  {
    _U[elem_j][0] = _U[elem_j][0]/_L[elem_j];
    _U[elem_j][1] = _U[elem_j][1] - _Dt_on_Dx*(_Pi_star[elem_j+1] - _Pi_star[elem_j]);
    //rho*v reste constant lors de l'étape acoustique
    _U[elem_j][3] = _U[elem_j][3] - _Dt_on_Dx*(_Pi_star[elem_j+1]*_u_star[elem_j+1] - _Pi_star[elem_j]*_u_star[elem_j]);
  }

  std::cout << "fin de l'étape acoustique, u = " << std::endl;
  for (int elem_j = 0; elem_j < _nbr_elements; elem_j++)
  {
    std::cout << _U[elem_j][1]/_U[elem_j][0] << std::endl;
  }

}

void Probleme1D::SaveIteration(int time_it)
{}
