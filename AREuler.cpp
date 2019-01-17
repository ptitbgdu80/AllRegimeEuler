#include "AREuler.h"

//Constructeur
Probleme1D::Probleme1D(int nbr_elements, double t_final, std::string file_name)
{
  _nbr_elements = nbr_elements;
  _t_final = t_final;
  _file_name = file_name;
  _delta_x = 1.0/nbr_elements;

  _cfl = 0.1;

  _U.resize(nbr_elements); // U contient les variables conservatives (rho, rho*u, rho*v, rho*E)
  _u.resize(nbr_elements);
  _a.resize(nbr_elements + 1);
  _u_star.resize(nbr_elements + 1); //_u_star[j] = u*_(j-1/2) (u* sur l'interface gauche de l'élément contenant xj)
  _L.resize(nbr_elements); // L[j] = 1 + Dt/Dx*(u*_(j+1/2) - u*_(j-1/2))
  _Pi.resize(nbr_elements); // Pi contient les pressions calculées à partir de U
  _Pi_star.resize(nbr_elements + 1); //_Pi_star[j] = Pi*_(j-1/2) (Pi* sur l'interface gauche de l'élément contenant xj)
  _Phi_interface.resize(nbr_elements + 1);


  //Conditions initiales
  for (int elem_j = 0; elem_j < nbr_elements; elem_j++)
  {
    double rho, u, v, p;
    if ((elem_j + 1./2.)*_delta_x < 0.5)
    {
      rho = 1.0;
      u = 0.0;
      v = 0.0;
      p = 1.0e5;
      _Pi[elem_j] = p;
    }

    else
    {
      rho = 0.1;
      u = 0.0;
      v = 0.0;
      p = 1.0e4;
      _Pi[elem_j] = p;
    }

    _U[elem_j].resize(4);
    _Phi_interface[elem_j].resize(4);

    _U[elem_j][0] = rho;
    _U[elem_j][1] = rho*u;
    _U[elem_j][2] = rho*v;
    _U[elem_j][3] = PressureToRhoE(p, _U[elem_j][0], _U[elem_j][1], _U[elem_j][2]);
  }


  //Conditions aux bords
  _left_bound_U = LeftBoundValue();
  _right_bound_U = RightBoundValue();

  _left_bound_u = _left_bound_U[1]/_left_bound_U[0];
  _right_bound_u = _right_bound_U[1]/_right_bound_U[0];
}



//Destructeur
Probleme1D::~Probleme1D()
{}



//Methodes de la classe
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
    _L[elem_j] = 1.0 + _Dt_on_Dx*(_u_star[elem_j+1] - _u_star[elem_j]);
  }
}

void Probleme1D::Update_Phi_interface()
{
  if (_u_star[0] >= 0)
  {
    _Phi_interface[0] = _left_bound_U;
  }
  else
  {
    _Phi_interface[0] = _U[0];
  }

  for (int elem_j = 1; elem_j < _nbr_elements; elem_j++)
  {
    if (_u_star[elem_j] >= 0)
    {
      _Phi_interface[elem_j] = _U[elem_j-1];
    }
    else
    {
      _Phi_interface[elem_j] = _U[elem_j];
    }
  }

  if (_u_star[_nbr_elements] >= 0)
  {
    _Phi_interface[_nbr_elements] = _U[_nbr_elements-1];
  }
  else
  {
    _Phi_interface[0] = _right_bound_U;
  }
}

std::vector<double> Probleme1D::LeftBoundValue()
{
  double rho, rho_u, rho_v, rho_E;
  // rho = 1.0;
  // u = 5000.0;
  // v = 5.0;
  //
  // rho_u = rho*u;
  // rho_v = rho*v;
  // rho_E = PressureToRhoE(_left_bound_Pi, rho, rho_u, rho_v);

  _left_bound_Pi = _Pi[0];

  rho = _U[0][0];
  rho_u = _U[0][1];
  rho_v = _U[0][2];
  rho_E = PressureToRhoE(_left_bound_Pi, rho, rho_u, rho_v);
  return {rho, rho_u, rho_v, rho_E};
}

std::vector<double> Probleme1D::RightBoundValue()
{
  double rho, rho_u, rho_v, rho_E;
  // rho = 1.0;
  // u = 5000.0;
  // v = 5.0;
  //
  // rho_u = rho*u;
  // rho_v = rho*v;
  // rho_E = PressureToRhoE(_right_bound_Pi, rho, rho_u, rho_v);

  _right_bound_Pi = _Pi[_nbr_elements-1];

  rho = _U[_nbr_elements-1][0];
  rho_u = _U[_nbr_elements-1][1];
  rho_v = _U[_nbr_elements-1][2];
  rho_E = PressureToRhoE(_right_bound_Pi, rho, rho_u, rho_v);
  return {rho, rho_u, rho_v, rho_E};
}

void Probleme1D::AcousticStep()
{
  for (int elem_j = 0; elem_j < _nbr_elements; elem_j++)
  {
    _U[elem_j][0] = _U[elem_j][0]/_L[elem_j];
    _U[elem_j][1] = (_U[elem_j][1] - _Dt_on_Dx*(_Pi_star[elem_j+1] - _Pi_star[elem_j]))/_L[elem_j];
    _U[elem_j][2] = _U[elem_j][2]/_L[elem_j];
    _U[elem_j][3] = (_U[elem_j][3] - _Dt_on_Dx*(_Pi_star[elem_j+1]*_u_star[elem_j+1] - _Pi_star[elem_j]*_u_star[elem_j]))/_L[elem_j];
  }
}

void Probleme1D::TransportStep()
{
  for (int elem_j = 0; elem_j < _nbr_elements; elem_j++)
  {
    _U[elem_j][0] = _U[elem_j][0]*_L[elem_j] + _Dt_on_Dx*(_u_star[elem_j+1]*_Phi_interface[elem_j+1][0] - _u_star[elem_j]*_Phi_interface[elem_j][0]);
    _U[elem_j][1] = _U[elem_j][1]*_L[elem_j] + _Dt_on_Dx*(_u_star[elem_j+1]*_Phi_interface[elem_j+1][1] - _u_star[elem_j]*_Phi_interface[elem_j][1]);
    _U[elem_j][2] = _U[elem_j][2]*_L[elem_j] + _Dt_on_Dx*(_u_star[elem_j+1]*_Phi_interface[elem_j+1][2] - _u_star[elem_j]*_Phi_interface[elem_j][2]);
    _U[elem_j][3] = _U[elem_j][3]*_L[elem_j] + _Dt_on_Dx*(_u_star[elem_j+1]*_Phi_interface[elem_j+1][3] - _u_star[elem_j]*_Phi_interface[elem_j][3]);
  }
}

void Probleme1D::SaveIteration(int time_it)
{
  if (time_it == 0)
  {
    system(("rm -rf " + _file_name).c_str());
    system(("mkdir " + _file_name).c_str());
  }

  std::ofstream file;
  file.open(_file_name + "/" + _file_name + std::to_string(time_it), std::ios::out);

  file << "# état du système à t = " << _time << std::endl;
  file << "-1";
  for (int iVar = 0; iVar < 4; iVar++)
  {
    file << " " << _left_bound_U[iVar];
  }
  file << " " << _left_bound_Pi << std::endl;

  for (int elem_j = 0; elem_j < _nbr_elements; elem_j++)
  {
    file << elem_j;
    for (int iVar = 0; iVar < 4; iVar++)
    {
      file << " " << _U[elem_j][iVar];
    }
    file << " " << _Pi[elem_j] << std::endl;
  }

  file << _nbr_elements;
  for (int iVar = 0; iVar < 4; iVar++)
  {
    file << " " << _right_bound_U[iVar];
  }
  file << " " << _right_bound_Pi << std::endl;

  // file << " " << _left_bound_U[0];
  // for (int iVar = 1; iVar < 4; iVar++)
  // {
  //   file << " " << _left_bound_U[iVar]/_left_bound_U[0];
  // }
  // file << " " << _left_bound_Pi << std::endl;
  //
  // for (int elem_j = 0; elem_j < _nbr_elements; elem_j++)
  // {
  //   file << elem_j;
  //   file << " " << _U[elem_j][0];
  //   for (int iVar = 1; iVar < 4; iVar++)
  //   {
  //     file << " " << _U[elem_j][iVar]/_U[elem_j][0];
  //   }
  //   file << " " << _Pi[elem_j] << std::endl;
  // }
  //
  // file << _nbr_elements;
  // file << " " << _right_bound_U[0];
  // for (int iVar = 1; iVar < 4; iVar++)
  // {
  //   file << " " << _right_bound_U[iVar]/_right_bound_U[0];
  // }
  // file << " " << _right_bound_Pi << std::endl;


  file.close();
}

void Probleme1D::TimeIteration(int time_it)
{
  Update_a();
  Update_u_star();
  Update_Pi_star();

  //Evaluation de la cfl pour l'acoustic step
  double max_for_cfl = _a[0]/std::min(_left_bound_U[0],_U[0][0]);
  double buffer;

  for (int elem_j = 1; elem_j < _nbr_elements; elem_j++)
  {
    buffer = _a[elem_j]/std::min(_U[elem_j-1][0],_U[elem_j][0]);
    if (buffer > max_for_cfl)
    {
      max_for_cfl = buffer;
    }
  }
  buffer = _a[_nbr_elements]/std::min(_U[_nbr_elements-1][0],_right_bound_U[0]);
  if (buffer > max_for_cfl)
  {
    max_for_cfl = buffer;
  }

  _Dt_on_Dx = _cfl/(2*max_for_cfl);


  //Evaluation de la cfl pour la transport step
  max_for_cfl = 0.5*(_u_star[0] + abs(_u_star[0]) - _u_star[1] + abs(_u_star[1]));

  for (int elem_j = 1; elem_j < _nbr_elements; elem_j++)
  {
    double buffer = 0.5*(_u_star[elem_j] + abs(_u_star[elem_j]) - _u_star[elem_j+1] + abs(_u_star[elem_j+1]));
    if (buffer > max_for_cfl)
    {
      max_for_cfl = buffer;
    }
  }

  _Dt_on_Dx = std::min(_cfl/max_for_cfl,_Dt_on_Dx);

  _time += _Dt_on_Dx*_delta_x;

  Update_L();

  // std::cout << "u500 : " << _u[500] << std::endl;
  // std::cout << "pi500 : " << _Pi[500] << std::endl;
  // std::cout << "a500 : " << _a[500] << std::endl;
  // std::cout << "a501 : " << _a[501] << std::endl;
  // std::cout << "pi*500 : " << _Pi_star[500] << std::endl;
  // std::cout << "pi*501 : " << _Pi_star[501] << std::endl;
  // std::cout << "u*500 : " << _u_star[500] << std::endl;
  // std::cout << "u*501 : " << _u_star[501] << std::endl;
  // std::cout << "L500 : " << _L[500] << std::endl;

  AcousticStep();

  // std::cout << "rho500 : " << _U[500][0] << std::endl;
  // std::cout << "rho_u500 : " << _U[500][1] << std::endl;
  // std::cout << "rho_E500 : " << _U[500][3] << std::endl;

  Update_Phi_interface();

  TransportStep();

  Update_u();
  Update_Pi();

  // std::cout << "p500 : " << _Pi[500] << std::endl;

  SaveIteration(time_it);

  _left_bound_U = LeftBoundValue();
  _right_bound_U = RightBoundValue();
}

void Probleme1D::Solve()
{
  _time = 0.0;
  int time_it = 0;

  Update_u();
  Update_Pi();

  SaveIteration(time_it);

  while (_time < _t_final)
  {
    std::cout << "t = " << _time << std::endl;
    time_it += 1;
    TimeIteration(time_it);
  }
}


//Méthodes indépendantes
double PressureToRhoE(double p, double rho, double rho_u, double rho_v)
{
  double rho_E = 2.5*p + 0.5*(rho_u*rho_u + rho_v*rho_v)/rho;
  return rho_E;
}
