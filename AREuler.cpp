#include "AREuler.h"

//Constructeur
Problem1D::Problem1D(int nbr_elements, double t_final, int choice_theta, std::string file_name)
{
  _nbr_elements = nbr_elements;
  _t_final = t_final;
  _choice_theta = choice_theta;
  _file_name = file_name;
  _delta_x = 1.0/nbr_elements;

  _cfl = 0.95;

  _U.resize(nbr_elements); // U contient les variables conservatives (rho, rho*u, rho*v, rho*E)
  _u.resize(nbr_elements);
  _a.resize(nbr_elements + 1);
  _theta.resize(nbr_elements + 1);
  _u_star.resize(nbr_elements + 1); //_u_star[j] = u*_(j-1/2) (u* sur l'interface gauche de l'élément contenant xj)
  _L.resize(nbr_elements); // L[j] = 1 + Dt/Dx*(u*_(j+1/2) - u*_(j-1/2))
  _Pi.resize(nbr_elements); // Pi contient les pressions calculées à partir de U
  _Pi_star.resize(nbr_elements + 1); //_Pi_star[j] = Pi*_(j-1/2) (Pi* sur l'interface gauche de l'élément contenant xj)
  _Phi_interface.resize(nbr_elements + 1);
  _left_bound_U.resize(4);
  _right_bound_U.resize(4);

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
    }

    else
    {
      rho = 0.125;
      u = 0.0;
      v = 0.0;
      p = 1.0e4;
    }

    _u[elem_j] = u;
    _Pi[elem_j] = p;

    _U[elem_j].resize(4);
    _Phi_interface[elem_j].resize(4);

    _U[elem_j][0] = rho;
    _U[elem_j][1] = rho*u;
    _U[elem_j][2] = rho*v;
    _U[elem_j][3] = PressureToRhoE(p, _U[elem_j][0], _U[elem_j][1], _U[elem_j][2]);
  }


  //Conditions aux bords
  Update_CL();

  //Initialisation des theta
  switch (choice_theta)
  {
    case 0 :
    for (int elem_j = 0; elem_j < nbr_elements+1; elem_j++)
    {
      _theta[elem_j] = 0;
    }
    break;

    case 1 :
    for (int elem_j = 0; elem_j < nbr_elements+1; elem_j++)
    {
      _theta[elem_j] = 1;
    }
    break;

    case 2 :
    break;

    default :
    std::cout << "Le choix de theta doit être 0, 1 ou 2" << std::endl;
    exit(1);
  }
}



//Destructeur
Problem1D::~Problem1D()
{}



//Methodes de la classe
void Problem1D::Update_u()
{
  for (int elem_j = 0; elem_j < _nbr_elements; elem_j++)
  {
    _u[elem_j] = _U[elem_j][1]/_U[elem_j][0];
  }
}

void Problem1D::Update_Pi()
{
  for (int elem_j = 0; elem_j < _nbr_elements; elem_j++)
  {
    _Pi[elem_j] = 0.4*(_U[elem_j][3] - (_U[elem_j][1]*_U[elem_j][1] + _U[elem_j][2]*_U[elem_j][2])/(2*_U[elem_j][0])); // p = (gamma - 1)*rho*e avec gamma = 1.4 et e = E - |u²|
  }
}

void Problem1D::Update_a()
{
  // a = 1.5*max(rho_L*c_L, rho_R*c_R) , c = sqrt(gamma*p/rho) => rho*c = sqrt(gamma*p*rho)

  double K = 1.5;
  double rhoLcL = sqrt(1.4*_left_bound_Pi*_left_bound_U[0]);
  double rhoRcR = sqrt(1.4*_Pi[0]*_U[0][0]);
  _a[0] = K*std::max(rhoLcL, rhoRcR);

  for (int elem_j = 1; elem_j < _nbr_elements; elem_j++)
  {
    rhoLcL = rhoRcR;
    rhoRcR = sqrt(1.4*_Pi[elem_j]*_U[elem_j][0]);
    _a[elem_j] = K*std::max(rhoLcL, rhoRcR);
  }

  rhoLcL = rhoRcR;
  rhoRcR = sqrt(1.4*_right_bound_Pi*_right_bound_U[0]);
  _a[_nbr_elements] = K*std::max(rhoLcL, rhoRcR);
}

void Problem1D::Update_theta()
{
  double cL, cR;
  switch (_choice_theta)
  {
    case 0 :
    case 1 :
    break;

    case 2 :
    cL = sqrt(1.4*_left_bound_Pi/_left_bound_U[0]);
    cR = sqrt(1.4*_Pi[0]/_U[0][0]);
    _theta[0] = std::min(abs(_u_star[0])/std::max(cL,cR),1.0);

    for (int elem_j = 1; elem_j < _nbr_elements; elem_j++)
    {
      cL = cR;
      cR = sqrt(1.4*_Pi[elem_j]/_U[elem_j][0]);
      _theta[elem_j] = std::min(abs(_u_star[elem_j])/std::max(cL,cR),1.0);
    }

    cL = cR;
    cR = sqrt(1.4*_right_bound_Pi/_right_bound_U[0]);
    _theta[_nbr_elements] = std::min(abs(_u_star[_nbr_elements])/std::max(cL,cR),1.0);
    break;

    default:
    std::cout << "Le choix de theta doit être 0, 1 ou 2" << std::endl;
    exit(1);
  }
}

void Problem1D::Update_u_star()
{
  _u_star[0] = 0.5*(_u[0] + _left_bound_u - (_Pi[0] - _left_bound_Pi)/_a[0]);
  for (int elem_j = 1; elem_j < _nbr_elements; elem_j++)
  {
    _u_star[elem_j] = 0.5*(_u[elem_j] + _u[elem_j-1] - (_Pi[elem_j] - _Pi[elem_j-1])/_a[elem_j]);
  }
  _u_star[_nbr_elements] = 0.5*(_right_bound_u + _u[_nbr_elements-1] - (_right_bound_Pi - _Pi[_nbr_elements-1])/_a[_nbr_elements]);
}

void Problem1D::Update_Pi_star()
{
  _Pi_star[0] = 0.5*(_Pi[0] + _left_bound_Pi - _theta[0]*_a[0]*(_u[0] - _left_bound_u));
  for (int elem_j = 1; elem_j < _nbr_elements; elem_j++)
  {
    _Pi_star[elem_j] = 0.5*(_Pi[elem_j] + _Pi[elem_j-1] - _theta[elem_j]*_a[elem_j]*(_u[elem_j] - _u[elem_j-1]));
  }
  _Pi_star[_nbr_elements] = 0.5*(_right_bound_Pi + _Pi[_nbr_elements-1] - _theta[_nbr_elements]*_a[_nbr_elements]*(_right_bound_u - _u[_nbr_elements-1]));
}

void Problem1D::Update_L()
{
  for (int elem_j = 0; elem_j < _nbr_elements; elem_j++)
  {
    _L[elem_j] = 1.0 + _Dt_on_Dx*(_u_star[elem_j+1] - _u_star[elem_j]);
  }
}

void Problem1D::Update_Phi_interface()
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

void Problem1D::Update_CL()
{
  _left_bound_Pi = _Pi[0];
  _right_bound_Pi = _Pi[_nbr_elements-1];

  _left_bound_U[0] = _U[0][0];
  _left_bound_U[1] = _U[0][1];
  _left_bound_U[2] = _U[0][2];
  _left_bound_U[3] = PressureToRhoE(_left_bound_Pi,  _U[0][0],  _U[0][1],  _U[0][2]);

  _right_bound_U[0] = _U[_nbr_elements-1][0];
  _right_bound_U[1] = _U[_nbr_elements-1][1];
  _right_bound_U[2] = _U[_nbr_elements-1][2];
  _right_bound_U[3] = PressureToRhoE(_left_bound_Pi, _U[_nbr_elements-1][0], _U[_nbr_elements-1][1], _U[_nbr_elements-1][3]);

  _left_bound_u = _left_bound_U[1]/_left_bound_U[0];
  _right_bound_u = _right_bound_U[1]/_right_bound_U[0];
}

void Problem1D::AcousticStep()
{
  for (int elem_j = 0; elem_j < _nbr_elements; elem_j++)
  {
    _U[elem_j][0] = _U[elem_j][0]/_L[elem_j];
    _U[elem_j][1] = (_U[elem_j][1] - _Dt_on_Dx*(_Pi_star[elem_j+1] - _Pi_star[elem_j]))/_L[elem_j];
    _U[elem_j][2] = _U[elem_j][2]/_L[elem_j];
    _U[elem_j][3] = (_U[elem_j][3] - _Dt_on_Dx*(_Pi_star[elem_j+1]*_u_star[elem_j+1] - _Pi_star[elem_j]*_u_star[elem_j]))/_L[elem_j];
  }
}

void Problem1D::TransportStep()
{
  for (int elem_j = 0; elem_j < _nbr_elements; elem_j++)
  {
    _U[elem_j][0] = _U[elem_j][0]*_L[elem_j] - _Dt_on_Dx*(_u_star[elem_j+1]*_Phi_interface[elem_j+1][0] - _u_star[elem_j]*_Phi_interface[elem_j][0]);
    _U[elem_j][1] = _U[elem_j][1]*_L[elem_j] - _Dt_on_Dx*(_u_star[elem_j+1]*_Phi_interface[elem_j+1][1] - _u_star[elem_j]*_Phi_interface[elem_j][1]);
    _U[elem_j][2] = _U[elem_j][2]*_L[elem_j] - _Dt_on_Dx*(_u_star[elem_j+1]*_Phi_interface[elem_j+1][2] - _u_star[elem_j]*_Phi_interface[elem_j][2]);
    _U[elem_j][3] = _U[elem_j][3]*_L[elem_j] - _Dt_on_Dx*(_u_star[elem_j+1]*_Phi_interface[elem_j+1][3] - _u_star[elem_j]*_Phi_interface[elem_j][3]);
  }
}

void Problem1D::SaveIteration(int time_it)
{
  if (time_it == 0)
  {
    system(("rm -rf " + _file_name).c_str());
    system(("mkdir " + _file_name).c_str());
  }

  std::ofstream file;
  file.open(_file_name + "/it_" + std::to_string(time_it), std::ios::out);

  file << "# état du système à t = " << _time << " : x, densité, pression, vitesse scalaire, nombre de Mach" << std::endl;
  file << -0.5*_delta_x;
  file << " " << _left_bound_U[0];
  file << " " << _left_bound_Pi;
  double velocity_mag = sqrt(_left_bound_U[1]*_left_bound_U[1] + _left_bound_U[2]*_left_bound_U[2])/_left_bound_U[0];
  double sound_speed = sqrt(1.4*_left_bound_Pi/_left_bound_U[0]);
  file << " " << velocity_mag;
  file << " " << velocity_mag/sound_speed;

  file << " " << 0 << " " << 0 << " " << 0 << " " << 1;

  file << std::endl;

  for (int elem_j = 0; elem_j < _nbr_elements; elem_j++)
  {
    file << (elem_j + 0.5)*_delta_x;
    file << " " << _U[elem_j][0];
    file << " " << _Pi[elem_j];
    velocity_mag = sqrt(_U[elem_j][1]*_U[elem_j][1] + _U[elem_j][2]*_U[elem_j][2])/_U[elem_j][0];
    sound_speed = sqrt(1.4*_Pi[elem_j]/_U[elem_j][0]);
    file << " " << velocity_mag;
    file << " " << velocity_mag/sound_speed;

    file << " " << _U[elem_j][3];
    // file << " " << _a[elem_j] << " " << _u_star[elem_j] << " " << _Pi_star[elem_j] << " " << _L[elem_j];
    file << std::endl;
  }

  file << 1.0 + 0.5*_delta_x;
  file << " " << _right_bound_U[0];
  file << " " << _right_bound_Pi;
  velocity_mag = sqrt(_right_bound_U[1]*_right_bound_U[1] + _right_bound_U[2]*_right_bound_U[2])/_right_bound_U[0];
  sound_speed = sqrt(1.4*_right_bound_Pi/_right_bound_U[0]);
  file << " " << velocity_mag;
  file << " " << velocity_mag/sound_speed;

  file << " " << _a[_nbr_elements] << " " << _u_star[_nbr_elements] << " " << _Pi_star[_nbr_elements] << " " << 1;

  file << std::endl;

  file.close();
}

void Problem1D::TimeIteration(int time_it)
{
  Update_a();
  Update_u_star();
  Update_theta();
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
    buffer = 0.5*(_u_star[elem_j] + abs(_u_star[elem_j]) - _u_star[elem_j+1] + abs(_u_star[elem_j+1]));
    if (buffer > max_for_cfl)
    {
      max_for_cfl = buffer;
    }
  }

  _Dt_on_Dx = std::min(_cfl/max_for_cfl,_Dt_on_Dx);

  _time += _Dt_on_Dx*_delta_x;

  Update_L();

  AcousticStep();

  Update_Phi_interface();

  TransportStep();

  Update_u();
  Update_Pi();
  Update_CL();

  SaveIteration(time_it);
}

void Problem1D::Solve()
{
  _time = 0.0;
  int time_it = 0;

  SaveIteration(time_it);

  while (_time < _t_final)
  {
    time_it += 1;
    TimeIteration(time_it);
  }

  std::cout << "t_final = " << _t_final << " atteint en " << time_it << " itérations" << std::endl;
}


//Méthodes indépendantes
double PressureToRhoE(double p, double rho, double rho_u, double rho_v)
{
  double rho_E = 2.5*p + 0.5*(rho_u*rho_u + rho_v*rho_v)/rho;
  return rho_E;
}
