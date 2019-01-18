#include "AREuler.h"

//Constructeur
Probleme2D::Probleme2D(int nbr_elements_1D, double t_final, int choix_theta, std::string file_name)
{
  _nbr_elements_1D = nbr_elements_1D;
  _t_final = t_final;
  _choix_theta = choix_theta;
  _file_name = file_name;
  _delta_s = 1.0/nbr_elements_1D;

  _cfl = 0.95;

  _U.resize(nbr_elements_1D); // U contient les variables conservatives (rho, rho*u, rho*v, rho*E)
  for (int line = 0; line < nbr_elements_1D; line++)
  {
    _U[line].resize(nbr_elements_1D);
    for (int column = 0; column < nbr_elements_1D; column++)
    {
      _U[line][column].resize(4);
    }
  }

  _u.resize(_nbr_elements_1D);
  _v.resize(_nbr_elements_1D);
  _Pi.resize(_nbr_elements_1D); // Pi contient les pressions calculées à partir de U
  _L.resize(_nbr_elements_1D); // L[j] = 1 + Dt/Ds*somme(u*jk)
  for (int line = 0; line < nbr_elements_1D; line++)
  {
    _u[line].resize(nbr_elements_1D);
    _v[line].resize(nbr_elements_1D);
    _Pi[line].resize(nbr_elements_1D);
    _L[line].resize(nbr_elements_1D);
  }

  _a_LR.resize(_nbr_elements_1D);
  _theta_LR.resize(_nbr_elements_1D);
  _u_star.resize(_nbr_elements_1D);
  _Pi_star_LR.resize(_nbr_elements_1D);
  _Phi_interface_LR.resize(_nbr_elements_1D);
  for (int line = 0; line < nbr_elements_1D; line++)
  {
    _a_LR[line].resize(_nbr_elements_1D+1);
    _theta_LR[line].resize(_nbr_elements_1D+1);
    _u_star[line].resize(_nbr_elements_1D+1);
    _Pi_star_LR[line].resize(_nbr_elements_1D+1);
    _Phi_interface_LR[line].resize(_nbr_elements_1D+1);
    for (int column = 0; column < nbr_elements_1D+1; column++)
    {
      _Phi_interface_LR[line][column].resize(4);
    }
  }

  _a_DU.resize(_nbr_elements_1D+1);
  _theta_DU.resize(_nbr_elements_1D+1);
  _v_star.resize(_nbr_elements_1D+1);
  _Pi_star_DU.resize(_nbr_elements_1D+1);
  _Phi_interface_DU.resize(_nbr_elements_1D+1);
  for (int line = 0; line < nbr_elements_1D+1; line++)
  {
    _a_DU[line].resize(_nbr_elements_1D);
    _theta_DU[line].resize(_nbr_elements_1D);
    _v_star[line].resize(_nbr_elements_1D);
    _Pi_star_DU[line].resize(_nbr_elements_1D);
    _Phi_interface_DU[line].resize(_nbr_elements_1D);
    for (int column = 0; column < nbr_elements_1D; column++)
    {
      _Phi_interface_DU[line][column].resize(4);
    }
  }


  //Conditions initiales
  double pi = 3.141592653589793;
  for (int line = 0; line < nbr_elements_1D; line++)
  {
    for (int column = 0; column < nbr_elements_1D; column++)
    {
      double x ,y ,rho, u, v, p;
      x = (column + 0.5)*_delta_s;
      y = (line + 0.5)*_delta_s;
      rho = 1.0 - 0.5*tanh(y - 0.5);
      u = 2*pow(sin(pi*x),2)*sin(pi*y)*cos(pi*y);
      v = -2*sin(pi*x)*cos(pi*x)*pow(sin(pi*y),2);
      p = 1.0e3;
      _Pi[line][column] = p;

      _U[line][column][0] = rho;
      _U[line][column][1] = rho*u;
      _U[line][column][2] = rho*v;
      _U[line][column][3] = PressureToRhoE2D(p, _U[line][column][0], _U[line][column][1], _U[line][column][2]);
    }
  }


  //Conditions aux bords
  _left_bound_U = LeftBoundValue();
  _right_bound_U = RightBoundValue();
  _down_bound_U = DownBoundValue();
  _up_bound_U = UpBoundValue();

  for (int index = 0; index < _nbr_elements_1D; index++)
  {
    _left_bound_u[index] = _left_bound_U[index][1]/_left_bound_U[index][0];
    _right_bound_u[index] = _right_bound_U[index][1]/_right_bound_U[index][0];

    _down_bound_v[index] = _down_bound_U[index][2]/_down_bound_U[index][0];
    _up_bound_v[index] = _up_bound_U[index][2]/_up_bound_U[index][0];
  }


  //Initialisation des theta
  switch (choix_theta)
  {
    case 0 :
    for (int line = 0; line < nbr_elements_1D; line++)
    {
      for (int column = 0; column < nbr_elements_1D+1; column++)
      {
        _theta_LR[line][column] = 0.0;
      }
    }

    for (int column = 0; column < nbr_elements_1D; column++)
    {
      for (int line = 0; line < nbr_elements_1D+1; line++)
      {
        _theta_DU[line][column] = 0.0;
      }
    }
    break;

    case 1 :
    for (int line = 0; line < nbr_elements_1D; line++)
    {
      for (int column = 0; column < nbr_elements_1D+1; column++)
      {
        _theta_LR[line][column] = 1.0;
      }
    }

    for (int column = 0; column < nbr_elements_1D; column++)
    {
      for (int line = 0; line < nbr_elements_1D+1; line++)
      {
        _theta_DU[line][column] = 1.0;
      }
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
Probleme2D::~Probleme2D()
{}



//Methodes de la classe
void Probleme2D::Update_u_v_Pi()
{
  for (int line = 0; line < _nbr_elements_1D; line++)
  {
    for (int column = 0; column < _nbr_elements_1D; column++)
    {
      _u[line][column] = _U[line][column][1]/_U[line][column][0];
      _v[line][column] = _U[line][column][2]/_U[line][column][0];
      _Pi[line][column] = 0.4*(_U[line][column][3] - 0.5*(_u[line][column]*_u[line][column] + _v[line][column]*_v[line][column]));
      if (_Pi[line][column] < 0)
      {
        std::cout << "Attention pression négative pour l'élément " << line << ", " << column << std::endl;
      }
    }
  }
}

void Probleme2D::Update_a_u_v_star()
{
  for (int line = 0; line < _nbr_elements_1D; line++)
  {
    double rhoLcL = sqrt(1.4*_left_bound_Pi[line]*_left_bound_U[line][0]);
    double rhoRcR = sqrt(1.4*_Pi[line][0]*_U[line][0][0]);
    _a_LR[line][0] = 1.5*std::max(rhoLcL, rhoRcR);
    _u_star[line][0] = 0.5*(_u[line][0] + _left_bound_u[line] - (_Pi[line][0] - _left_bound_Pi[line])/_a_LR[line][0]);

    for (int column = 1; column < _nbr_elements_1D; column++)
    {
      rhoLcL = rhoRcR;
      rhoRcR = sqrt(1.4*_Pi[line][column]*_U[line][column][0]);
      _a_LR[line][column] = 1.5*std::max(rhoLcL, rhoRcR);
      _u_star[line][column] = 0.5*(_u[line][column] + _u[line][column-1] - (_Pi[line][column] - _Pi[line][column-1])/_a_LR[line][column]);
    }

    rhoLcL = rhoRcR;
    rhoRcR = sqrt(1.4*_right_bound_Pi[line]*_right_bound_U[line][0]);
    _a_LR[line][_nbr_elements_1D] = 1.5*std::max(rhoLcL, rhoRcR);
    _u_star[line][_nbr_elements_1D] = 0.5*(_right_bound_u[line] + _u[line][_nbr_elements_1D-1] - (_right_bound_Pi[line] - _Pi[line][_nbr_elements_1D-1])/_a_LR[line][_nbr_elements_1D]);
  }



  for (int column = 0; column < _nbr_elements_1D; column++)
  {
    double rhoDcD = sqrt(1.4*_down_bound_Pi[column]*_down_bound_U[column][0]);
    double rhoUcU = sqrt(1.4*_Pi[0][column]*_U[0][column][0]);
    _a_DU[0][column] = 1.5*std::max(rhoDcD, rhoUcU);
    _v_star[0][column] = 0.5*(_v[0][column] + _down_bound_v[column] - (_Pi[0][column] - _down_bound_Pi[column])/_a_DU[0][column]);

    for (int line = 1; line < _nbr_elements_1D; line++)
    {
      rhoDcD = rhoUcU;
      rhoUcU = sqrt(1.4*_Pi[line][column]*_U[line][column][0]);
      _a_DU[line][column] = 1.5*std::max(rhoDcD, rhoUcU);
      _v_star[line][column] = 0.5*(_v[line][column] + _v[line-1][column] - (_Pi[line][column] - _Pi[line-1][column])/_a_DU[line][column]);
    }

    rhoDcD = rhoUcU;
    rhoUcU = sqrt(1.4*_up_bound_Pi[column]*_up_bound_U[column][0]);
    _a_DU[_nbr_elements_1D][column] = 1.5*std::max(rhoDcD, rhoUcU);
    _v_star[_nbr_elements_1D][column] = 0.5*(_up_bound_v[column] + _v[_nbr_elements_1D-1][column] - (_up_bound_Pi[column] - _Pi[_nbr_elements_1D-1][column])/_a_DU[_nbr_elements_1D][column]);
  }
}

void Probleme2D::Update_theta()
{
  double cL, cR, cD, cU;
  switch (_choix_theta)
  {
    case 0 :
    case 1 :
    break;

    case 2 :
    for (int line = 0; line < _nbr_elements_1D; line++)
    {
      cL = sqrt(1.4*_left_bound_Pi[line]/_left_bound_U[line][0]);
      cR = sqrt(1.4*_Pi[line][0]/_U[line][0][0]);
      _theta_LR[line][0] = std::min(abs(_u_star[line][0])/std::max(cL,cR),1.0);

      for (int column = 1; column < _nbr_elements_1D; column++)
      {
        cL = cR;
        cR = sqrt(1.4*_Pi[line][column]/_U[line][column][0]);
        _theta_LR[line][column] = std::min(abs(_u_star[line][column])/std::max(cL,cR),1.0);
      }

      cL = cR;
      cR = sqrt(1.4*_right_bound_Pi[line]/_right_bound_U[line][0]);
      _theta_LR[line][_nbr_elements_1D] = std::min(abs(_u_star[line][_nbr_elements_1D])/std::max(cL,cR),1.0);
    }


    for (int column = 0; column < _nbr_elements_1D; column++)
    {
      cD = sqrt(1.4*_down_bound_Pi[column]/_down_bound_U[column][0]);
      cU = sqrt(1.4*_Pi[0][column]/_U[0][column][0]);
      _theta_DU[0][column] = std::min(abs(_v_star[0][column])/std::max(cD,cU),1.0);

      for (int line = 1; line < _nbr_elements_1D; line++)
      {
        cD = cU;
        cU = sqrt(1.4*_Pi[line][column]/_U[line][column][0]);
        _theta_DU[line][column] = std::min(abs(_v_star[line][column])/std::max(cD,cU),1.0);
      }

      cD = cU;
      cU = sqrt(1.4*_up_bound_Pi[column]/_up_bound_U[column][0]);
      _theta_DU[_nbr_elements_1D][column] = std::min(abs(_v_star[_nbr_elements_1D][column])/std::max(cD,cU),1.0);
    }

    break;

    default:
    std::cout << "Le choix de theta doit être 0, 1 ou 2" << std::endl;
    exit(1);
  }
}

void Probleme2D::Update_Pi_star()
{
  for (int line = 0; line < _nbr_elements_1D; line++)
  {
    _Pi_star_LR[line][0] = 0.5*(_Pi[line][0] + _left_bound_Pi[line] - _theta_LR[line][0]*_a_LR[line][0]*(_u[line][0] - _left_bound_u[line]));
    for (int column = 1; column < _nbr_elements_1D; column++)
    {
      _Pi_star_LR[line][column] = 0.5*(_Pi[line][column] + _Pi[line][column-1] - _theta_LR[line][column]*_a_LR[line][column]*(_u[line][column] - _u[line][column-1]));
    }
    _Pi_star_LR[line][_nbr_elements_1D] = 0.5*(_right_bound_Pi[line] + _Pi[line][_nbr_elements_1D-1] - _theta_LR[line][_nbr_elements_1D-1]*_a_LR[line][_nbr_elements_1D-1]*(_right_bound_u[line] - _u[line][_nbr_elements_1D-1]));
  }

  for (int column = 0; column < _nbr_elements_1D; column++)
  {
    _Pi_star_DU[0][column] = 0.5*(_Pi[0][column] + _down_bound_Pi[column] - _theta_DU[0][column]*_a_DU[0][column]*(_v[0][column] - _down_bound_v[column]));
    for (int line = 1; line < _nbr_elements_1D; line++)
    {
      _Pi_star_DU[line][column] = 0.5*(_Pi[line][column] + _Pi[line-1][column] - _theta_DU[line][column]*_a_DU[line][column]*(_v[line][column] - _v[line-1][column]));
    }
    _Pi_star_DU[_nbr_elements_1D][column] = 0.5*(_up_bound_Pi[column] + _Pi[_nbr_elements_1D-1][column] - _theta_DU[_nbr_elements_1D-1][column]*_a_LR[_nbr_elements_1D-1][column]*(_up_bound_v[column] - _u[_nbr_elements_1D-1][column]));
  }
}

void Probleme2D::Update_L()
{
  for (int line = 0; line < _nbr_elements_1D; line++)
  {
    for (int column = 0; column < _nbr_elements_1D; column++)
    {
      _L[line][column] = 1.0 + _Dt_on_Ds*(_u_star[line][column+1] - _u_star[line][column] + _v_star[line+1][column] - _v_star[line][column]);
    }
  }
}

void Probleme2D::Update_Phi_interface()
{
  for (int line = 0; line < _nbr_elements_1D; line++)
  {
    if (_u_star[line][0] >= 0)
    {
      _Phi_interface_LR[line][0] = _left_bound_U[line];
    }
    else
    {
      _Phi_interface_LR[line][0] = _U[line][0];
    }

    for (int column = 1; column < _nbr_elements_1D; column++)
    {
      if (_u_star[line][column] >= 0)
      {
        _Phi_interface_LR[line][column] = _U[line][column-1];
      }
      else
      {
        _Phi_interface_LR[line][column] = _U[line][column];
      }
    }

    if (_u_star[line][_nbr_elements_1D] >= 0)
    {
      _Phi_interface_LR[line][_nbr_elements_1D] = _U[line][_nbr_elements_1D-1];
    }
    else
    {
      _Phi_interface_LR[line][_nbr_elements_1D] = _right_bound_U[line];
    }
  }


  for (int column = 0; column < _nbr_elements_1D; column++)
  {
    if (_v_star[0][column] >= 0)
    {
      _Phi_interface_DU[0][column] = _down_bound_U[column];
    }
    else
    {
      _Phi_interface_DU[0][column] = _U[0][column];
    }

    for (int line = 1; line < _nbr_elements_1D; line++)
    {
      if (_v_star[line][column] >= 0)
      {
        _Phi_interface_DU[line][column] = _U[line-1][column];
      }
      else
      {
        _Phi_interface_DU[line][column] = _U[line][column];
      }
    }

    if (_v_star[_nbr_elements_1D][column] >= 0)
    {
      _Phi_interface_DU[_nbr_elements_1D][column] = _U[_nbr_elements_1D-1][column];
    }
    else
    {
      _Phi_interface_DU[_nbr_elements_1D][column] = _up_bound_U[column];
    }
  }
}
//
// std::vector<std::vector<double>> Probleme2D::LeftBoundValue()
// {
//   double rho, rho_u, rho_v, rho_E;
//
//   _left_bound_Pi = _Pi[0];
//
//   rho = _U[0][0];
//   rho_u = _U[0][1];
//   rho_v = _U[0][2];
//   rho_E = PressureToRhoE(_left_bound_Pi, rho, rho_u, rho_v);
//   return {rho, rho_u, rho_v, rho_E};
// }
//
// std::vector<std::vector<double>> Probleme2D::RightBoundValue()
// {
//   double rho, rho_u, rho_v, rho_E;
//
//   _right_bound_Pi = _Pi[_nbr_elements-1];
//
//   rho = _U[_nbr_elements-1][0];
//   rho_u = _U[_nbr_elements-1][1];
//   rho_v = _U[_nbr_elements-1][2];
//   rho_E = PressureToRhoE(_right_bound_Pi, rho, rho_u, rho_v);
//   return {rho, rho_u, rho_v, rho_E};
// }

void Probleme2D::AcousticStep()
{
  for (int line = 0; line < _nbr_elements_1D; line++)
  {
    for (int column = 0; column < _nbr_elements_1D; column++)
    {
      _U[line][column][0] = _U[line][column][0]/_L[line][column];
      _U[line][column][1] = (_U[line][column][1] - _Dt_on_Ds*(_Pi_star_LR[line][column+1] - _Pi_star_LR[line][column]))/_L[line][column];
      _U[line][column][2] = (_U[line][column][2] - _Dt_on_Ds*(_Pi_star_DU[line-1][column] - _Pi_star_DU[line][column]))/_L[line][column];
      _U[line][column][3] = (_U[line][column][3] - _Dt_on_Ds*(_Pi_star_LR[line][column+1]*_u_star[line][column+1] - _Pi_star_LR[line][column]*_u_star[line][column] + _Pi_star_DU[line-1][column]*_v_star[line-1][column] - _Pi_star_DU[line][column]*_v_star[line][column]))/_L[line][column];
    }
  }
}
//
// void Probleme2D::TransportStep()
// {
//   for (int elem_j = 0; elem_j < _nbr_elements; elem_j++)
//   {
//     _U[elem_j][0] = _U[elem_j][0]*_L[elem_j] - _Dt_on_Dx*(_u_star[elem_j+1]*_Phi_interface[elem_j+1][0] - _u_star[elem_j]*_Phi_interface[elem_j][0]);
//     _U[elem_j][1] = _U[elem_j][1]*_L[elem_j] - _Dt_on_Dx*(_u_star[elem_j+1]*_Phi_interface[elem_j+1][1] - _u_star[elem_j]*_Phi_interface[elem_j][1]);
//     _U[elem_j][2] = _U[elem_j][2]*_L[elem_j] - _Dt_on_Dx*(_u_star[elem_j+1]*_Phi_interface[elem_j+1][2] - _u_star[elem_j]*_Phi_interface[elem_j][2]);
//     _U[elem_j][3] = _U[elem_j][3]*_L[elem_j] - _Dt_on_Dx*(_u_star[elem_j+1]*_Phi_interface[elem_j+1][3] - _u_star[elem_j]*_Phi_interface[elem_j][3]);
//   }
// }
//
// void Probleme2D::SaveIteration(int time_it)
// {
//   if (time_it == 0)
//   {
//     system(("rm -rf " + _file_name).c_str());
//     system(("mkdir " + _file_name).c_str());
//   }
//
//   std::ofstream file;
//   file.open(_file_name + "/" + _file_name + std::to_string(time_it), std::ios::out);
//
//   file << "# état du système à t = " << _time << std::endl;
//   file << -_delta_x/2.;
//   for (int iVar = 0; iVar < 4; iVar++)
//   {
//     file << " " << _left_bound_U[iVar];
//   }
//   file << " " << _left_bound_Pi << std::endl;
//
//   for (int elem_j = 0; elem_j < _nbr_elements; elem_j++)
//   {
//     file << (elem_j + 0.5)*_delta_x;
//     for (int iVar = 0; iVar < 4; iVar++)
//     {
//       file << " " << _U[elem_j][iVar];
//     }
//     file << " " << _Pi[elem_j] << std::endl;
//   }
//
//   file << 1.0 + 0.5*_delta_x;
//   for (int iVar = 0; iVar < 4; iVar++)
//   {
//     file << " " << _right_bound_U[iVar];
//   }
//   file << " " << _right_bound_Pi << std::endl;
//
//   // file << " " << _left_bound_U[0];
//   // for (int iVar = 1; iVar < 4; iVar++)
//   // {
//   //   file << " " << _left_bound_U[iVar]/_left_bound_U[0];
//   // }
//   // file << " " << _left_bound_Pi << std::endl;
//   //
//   // for (int elem_j = 0; elem_j < _nbr_elements; elem_j++)
//   // {
//   //   file << elem_j;
//   //   file << " " << _U[elem_j][0];
//   //   for (int iVar = 1; iVar < 4; iVar++)
//   //   {
//   //     file << " " << _U[elem_j][iVar]/_U[elem_j][0];
//   //   }
//   //   file << " " << _Pi[elem_j] << std::endl;
//   // }
//   //
//   // file << _nbr_elements;
//   // file << " " << _right_bound_U[0];
//   // for (int iVar = 1; iVar < 4; iVar++)
//   // {
//   //   file << " " << _right_bound_U[iVar]/_right_bound_U[0];
//   // }
//   // file << " " << _right_bound_Pi << std::endl;
//
//
//   file.close();
// }
//
// void Probleme2D::TimeIteration(int time_it)
// {
//   Update_a();
//   Update_u_star();
//   Update_Pi_star();
//   Update_theta();
//
//   //Evaluation de la cfl pour l'acoustic step
//   double max_for_cfl = _a[0]/std::min(_left_bound_U[0],_U[0][0]);
//   double buffer;
//
//   for (int elem_j = 1; elem_j < _nbr_elements; elem_j++)
//   {
//     buffer = _a[elem_j]/std::min(_U[elem_j-1][0],_U[elem_j][0]);
//     if (buffer > max_for_cfl)
//     {
//       max_for_cfl = buffer;
//     }
//   }
//   buffer = _a[_nbr_elements]/std::min(_U[_nbr_elements-1][0],_right_bound_U[0]);
//   if (buffer > max_for_cfl)
//   {
//     max_for_cfl = buffer;
//   }
//
//   _Dt_on_Dx = _cfl/(2*max_for_cfl);
//
//
//   //Evaluation de la cfl pour la transport step
//   max_for_cfl = 0.5*(_u_star[0] + abs(_u_star[0]) - _u_star[1] + abs(_u_star[1]));
//
//   for (int elem_j = 1; elem_j < _nbr_elements; elem_j++)
//   {
//     double buffer = 0.5*(_u_star[elem_j] + abs(_u_star[elem_j]) - _u_star[elem_j+1] + abs(_u_star[elem_j+1]));
//     if (buffer > max_for_cfl)
//     {
//       max_for_cfl = buffer;
//     }
//   }
//
//   _Dt_on_Dx = std::min(_cfl/max_for_cfl,_Dt_on_Dx);
//
//   _time += _Dt_on_Dx*_delta_x;
//
//   Update_L();
//
//   AcousticStep();
//
//   Update_Phi_interface();
//
//   TransportStep();
//
//   Update_u();
//   Update_Pi();
//
//   SaveIteration(time_it);
//
//   _left_bound_U = LeftBoundValue();
//   _right_bound_U = RightBoundValue();
// }
//
// void Probleme2D::Solve()
// {
//   _time = 0.0;
//   int time_it = 0;
//
//   Update_u();
//   Update_Pi();
//
//   SaveIteration(time_it);
//
//   while (_time < _t_final)
//   {
//     std::cout << "t = " << _time << std::endl;
//     time_it += 1;
//     TimeIteration(time_it);
//   }
// }


//Méthodes indépendantes
double PressureToRhoE2D(double p, double rho, double rho_u, double rho_v)
{
  double rho_E = 2.5*p + 0.5*(rho_u*rho_u + rho_v*rho_v)/rho;
  return rho_E;
}
