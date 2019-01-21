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

  _u.resize(_nbr_elements_1D);
  _v.resize(_nbr_elements_1D);
  _Pi.resize(_nbr_elements_1D); // Pi contient les pressions calculées à partir de U
  _L.resize(_nbr_elements_1D); // L[j] = 1 + Dt/Ds*somme(u*jk)
  _U.resize(nbr_elements_1D); // U contient les variables conservatives (rho, rho*u, rho*v, rho*E)
  _left_bound_Pi.resize(nbr_elements_1D);
  _right_bound_Pi.resize(nbr_elements_1D);
  _down_bound_Pi.resize(nbr_elements_1D);
  _up_bound_Pi.resize(nbr_elements_1D);
  _left_bound_u.resize(nbr_elements_1D);
  _right_bound_u.resize(nbr_elements_1D);
  _down_bound_v.resize(nbr_elements_1D);
  _up_bound_v.resize(nbr_elements_1D);
  _left_bound_U.resize(nbr_elements_1D);
  _right_bound_U.resize(nbr_elements_1D);
  _down_bound_U.resize(nbr_elements_1D);
  _up_bound_U.resize(nbr_elements_1D);

  for (int line = 0; line < nbr_elements_1D; line++)
  {
    _u[line].resize(nbr_elements_1D);
    _v[line].resize(nbr_elements_1D);
    _Pi[line].resize(nbr_elements_1D);
    _L[line].resize(nbr_elements_1D);
    _U[line].resize(nbr_elements_1D);
    _left_bound_U[line].resize(4);
    _right_bound_U[line].resize(4);
    _down_bound_U[line].resize(4); //C'est un parcours par colonne pour down et up mais flemme de les mettre ailleurs
    _up_bound_U[line].resize(4);

    for (int column = 0; column < nbr_elements_1D; column++)
    {
      _U[line][column].resize(4);
    }
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

  // double pi = 3.141592653589793;
  // for (int line = 0; line < nbr_elements_1D; line++)
  // {
  //   for (int column = 0; column < nbr_elements_1D; column++)
  //   {
  //     double x ,y ,rho, u, v, p;
  //     x = (column + 0.5)*_delta_s;
  //     y = (line + 0.5)*_delta_s;
  //     rho = 1.0 - 0.5*tanh(y - 0.5);
  //     u = 2*pow(sin(pi*x),2)*sin(pi*y)*cos(pi*y);
  //     v = -2*sin(pi*x)*cos(pi*x)*pow(sin(pi*y),2);
  //     p = 1.0e3;
  //     _u[line][column] = u;
  //     _v[line][column] = v;
  //     _Pi[line][column] = p;
  //
  //     _U[line][column][0] = rho;
  //     _U[line][column][1] = rho*u;
  //     _U[line][column][2] = rho*v;
  //     _U[line][column][3] = PressureToRhoE2D(p, _U[line][column][0], _U[line][column][1], _U[line][column][2]);
  //   }
  // }

  for (int line = 0; line < nbr_elements_1D; line++)
  {
    for (int column = 0; column < nbr_elements_1D; column++)
    {
      double x ,y ,rho, u, v, p;
      x = (column + 0.5)*_delta_s;
      y = (line + 0.5)*_delta_s;
      if (x < 0.5)
      {
        if (y < 0.5)
        {
          rho = 0.138;
          u = 1.206;
          v = 1.206;
          p = 0.029;
        }
        else
        {
          rho = 0.5323;
          u = 1.206;
          v = 0.0;
          p = 0.3;
        }
      }
      else if (y < 0.5)
      {
        rho = 0.5323;
        u = 0.0;
        v = 1.206;
        p = 0.3;
      }
      else
      {
        rho = 1.5;
        u = 0.0;
        v = 0.0;
        p = 1.5;
      }

      _u[line][column] = u;
      _v[line][column] = v;
      _Pi[line][column] = p;

      _U[line][column][0] = rho;
      _U[line][column][1] = rho*u;
      _U[line][column][2] = rho*v;
      _U[line][column][3] = PressureToRhoE2D(p, _U[line][column][0], _U[line][column][1], _U[line][column][2]);
    }
  }

  //Conditions aux bords
  Update_CL();

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
      _Pi[line][column] = 0.4*(_U[line][column][3] - 0.5*_U[line][column][0]*(_u[line][column]*_u[line][column] + _v[line][column]*_v[line][column]));
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
    _Pi_star_LR[line][_nbr_elements_1D] = 0.5*(_right_bound_Pi[line] + _Pi[line][_nbr_elements_1D-1] - _theta_LR[line][_nbr_elements_1D]*_a_LR[line][_nbr_elements_1D]*(_right_bound_u[line] - _u[line][_nbr_elements_1D-1]));
  }

  for (int column = 0; column < _nbr_elements_1D; column++)
  {
    _Pi_star_DU[0][column] = 0.5*(_Pi[0][column] + _down_bound_Pi[column] - _theta_DU[0][column]*_a_DU[0][column]*(_v[0][column] - _down_bound_v[column]));
    for (int line = 1; line < _nbr_elements_1D; line++)
    {
      _Pi_star_DU[line][column] = 0.5*(_Pi[line][column] + _Pi[line-1][column] - _theta_DU[line][column]*_a_DU[line][column]*(_v[line][column] - _v[line-1][column]));
    }
    _Pi_star_DU[_nbr_elements_1D][column] = 0.5*(_up_bound_Pi[column] + _Pi[_nbr_elements_1D-1][column] - _theta_DU[_nbr_elements_1D][column]*_a_DU[_nbr_elements_1D][column]*(_up_bound_v[column] - _v[_nbr_elements_1D-1][column]));
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

void Probleme2D::Update_CL()
{
  // for (int index = 0; index < _nbr_elements_1D; index++)
  // {
  //   _left_bound_Pi[index] = _Pi[index][0];
  //   _right_bound_Pi[index] = _Pi[index][_nbr_elements_1D-1];
  //   _down_bound_Pi[index] = _Pi[0][index];
  //   _up_bound_Pi[index] = _Pi[_nbr_elements_1D-1][index];
  //
  //   _left_bound_U[index][0] = _U[index][0][0];
  //   _left_bound_U[index][1] = 0.0;
  //   _left_bound_U[index][2] = 0.0;
  //   _left_bound_U[index][3] = PressureToRhoE(_left_bound_Pi[index], _left_bound_U[index][0], 0, 0);
  //
  //   _right_bound_U[index][0] = _U[index][_nbr_elements_1D-1][0];
  //   _right_bound_U[index][1] = 0.0;
  //   _right_bound_U[index][2] = 0.0;
  //   _right_bound_U[index][3] = PressureToRhoE(_right_bound_Pi[index], _right_bound_U[index][0], 0, 0);
  //
  //   _down_bound_U[index][0] = _U[0][index][0];
  //   _down_bound_U[index][1] = 0.0;
  //   _down_bound_U[index][2] = 0.0;
  //   _down_bound_U[index][3] = PressureToRhoE(_down_bound_Pi[index], _down_bound_U[index][0], 0, 0);
  //
  //   _up_bound_U[index][0] = _U[_nbr_elements_1D-1][index][0];
  //   _up_bound_U[index][1] = 0.0;
  //   _up_bound_U[index][2] = 0.0;
  //   _up_bound_U[index][3] = PressureToRhoE(_up_bound_Pi[index], _up_bound_U[index][0], 0, 0);
  //
  //   _left_bound_u[index] = 0.0;
  //   _right_bound_u[index] = 0.0;
  //   _down_bound_v[index] = 0.0;
  //   _up_bound_v[index] = 0.0;
  // }

  for (int index = 0; index < _nbr_elements_1D; index++)
  {
    _left_bound_Pi[index] = _Pi[index][0];
    _right_bound_Pi[index] = _Pi[index][_nbr_elements_1D-1];
    _down_bound_Pi[index] = _Pi[0][index];
    _up_bound_Pi[index] = _Pi[_nbr_elements_1D-1][index];

    _left_bound_U[index][0] = _U[index][0][0];
    _left_bound_U[index][1] = _U[index][0][1];
    _left_bound_U[index][2] = _U[index][0][2];
    _left_bound_U[index][3] = PressureToRhoE(_left_bound_Pi[index], _left_bound_U[index][0], _left_bound_U[index][1], _left_bound_U[index][2]);

    _right_bound_U[index][0] = _U[index][_nbr_elements_1D-1][0];
    _right_bound_U[index][1] = _U[index][_nbr_elements_1D-1][1];
    _right_bound_U[index][2] = _U[index][_nbr_elements_1D-1][2];
    _right_bound_U[index][3] = PressureToRhoE(_right_bound_Pi[index], _right_bound_U[index][0], _right_bound_U[index][1], _right_bound_U[index][2]);

    _down_bound_U[index][0] = _U[0][index][0];
    _down_bound_U[index][1] = _U[0][index][1];
    _down_bound_U[index][2] = _U[0][index][2];
    _down_bound_U[index][3] = PressureToRhoE(_down_bound_Pi[index], _down_bound_U[index][0], _down_bound_U[index][1], _down_bound_U[index][2]);

    _up_bound_U[index][0] = _U[_nbr_elements_1D-1][index][0];
    _up_bound_U[index][1] = _U[_nbr_elements_1D-1][index][1];
    _up_bound_U[index][2] = _U[_nbr_elements_1D-1][index][2];
    _up_bound_U[index][3] = PressureToRhoE(_up_bound_Pi[index], _up_bound_U[index][0], _up_bound_U[index][1], _up_bound_U[index][2]);

    _left_bound_u[index] = _u[index][0];
    _right_bound_u[index] = _u[index][_nbr_elements_1D-1];
    _down_bound_v[index] = _v[0][index];
    _up_bound_v[index] = _v[_nbr_elements_1D-1][index];
  }

}

void Probleme2D::AcousticStep()
{
  for (int line = 0; line < _nbr_elements_1D; line++)
  {
    for (int column = 0; column < _nbr_elements_1D; column++)
    {
      _U[line][column][0] = _U[line][column][0]/_L[line][column];
      _U[line][column][1] = (_U[line][column][1] - _Dt_on_Ds*(_Pi_star_LR[line][column+1] - _Pi_star_LR[line][column]))/_L[line][column];
      _U[line][column][2] = (_U[line][column][2] - _Dt_on_Ds*(_Pi_star_DU[line+1][column] - _Pi_star_DU[line][column]))/_L[line][column];
      _U[line][column][3] = (_U[line][column][3] - _Dt_on_Ds*(_Pi_star_LR[line][column+1]*_u_star[line][column+1] - _Pi_star_LR[line][column]*_u_star[line][column] + _Pi_star_DU[line+1][column]*_v_star[line+1][column] - _Pi_star_DU[line][column]*_v_star[line][column]))/_L[line][column];
    }
  }
}

void Probleme2D::TransportStep()
{
  for (int line = 0; line < _nbr_elements_1D; line++)
  {
    for (int column = 0; column < _nbr_elements_1D; column++)
    {
      _U[line][column][0] = _U[line][column][0]*_L[line][column];
      _U[line][column][0] -= _Dt_on_Ds*(_u_star[line][column+1]*_Phi_interface_LR[line][column+1][0] - _u_star[line][column]*_Phi_interface_LR[line][column][0]);
      _U[line][column][0] -= _Dt_on_Ds*(_v_star[line+1][column]*_Phi_interface_DU[line+1][column][0] - _v_star[line][column]*_Phi_interface_DU[line][column][0]);

      _U[line][column][1] = _U[line][column][1]*_L[line][column];
      _U[line][column][1] -= _Dt_on_Ds*(_u_star[line][column+1]*_Phi_interface_LR[line][column+1][1] - _u_star[line][column]*_Phi_interface_LR[line][column][1]);
      _U[line][column][1] -= _Dt_on_Ds*(_v_star[line+1][column]*_Phi_interface_DU[line+1][column][1] - _v_star[line][column]*_Phi_interface_DU[line][column][1]);

      _U[line][column][2] = _U[line][column][2]*_L[line][column];
      _U[line][column][2] -= _Dt_on_Ds*(_u_star[line][column+1]*_Phi_interface_LR[line][column+1][2] - _u_star[line][column]*_Phi_interface_LR[line][column][2]);
      _U[line][column][2] -= _Dt_on_Ds*(_v_star[line+1][column]*_Phi_interface_DU[line+1][column][2] - _v_star[line][column]*_Phi_interface_DU[line][column][2]);

      _U[line][column][3] = _U[line][column][3]*_L[line][column];
      _U[line][column][3] -= _Dt_on_Ds*(_u_star[line][column+1]*_Phi_interface_LR[line][column+1][3] - _u_star[line][column]*_Phi_interface_LR[line][column][3]);
      _U[line][column][3] -= _Dt_on_Ds*(_v_star[line+1][column]*_Phi_interface_DU[line+1][column][3] - _v_star[line][column]*_Phi_interface_DU[line][column][3]);
    }
  }
}

void Probleme2D::SaveIteration(int time_it)
{
  if (time_it == 0)
  {
    system(("rm -rf " + _file_name).c_str());
    system(("mkdir " + _file_name).c_str());
  }

  std::ofstream file;
  file.open(_file_name + "/" + _file_name + std::to_string(time_it), std::ios::out);

  file << "# état du système à t = " << _time << std::endl;
  file << -0.5*_delta_s << " " << -0.5*_delta_s;
  file << " " << _down_bound_U[0][0]; //Le coin a la même valeur que le premier élément du bas
  file << " " << _down_bound_Pi[0];
  double velocity_mag = sqrt(_down_bound_U[0][1]*_down_bound_U[0][1] + _down_bound_U[0][2]*_down_bound_U[0][2])/_down_bound_U[0][0];
  double sound_speed = sqrt(1.4*_down_bound_Pi[0]/_down_bound_U[0][0]);
  file << " " << velocity_mag;
  file << " " << velocity_mag/sound_speed;
  file << std::endl;

  for (int column = 0; column < _nbr_elements_1D; column++)
  {
    file << (column + 0.5)*_delta_s << " " << -0.5*_delta_s;
    file << " " << _down_bound_U[column][0];
    file << " " << _down_bound_Pi[column];
    velocity_mag = sqrt(_down_bound_U[column][1]*_down_bound_U[column][1] + _down_bound_U[column][2]*_down_bound_U[column][2])/_down_bound_U[column][0];
    sound_speed = sqrt(1.4*_down_bound_Pi[column]/_down_bound_U[column][0]);
    file << " " << velocity_mag;
    file << " " << velocity_mag/sound_speed;
    file << std::endl;
  }

  file << 1.0 + 0.5*_delta_s << " " << -0.5*_delta_s;
  file << " " << _down_bound_U[_nbr_elements_1D-1][0];
  file << " " << _down_bound_Pi[_nbr_elements_1D-1];
  velocity_mag = sqrt(_down_bound_U[_nbr_elements_1D-1][1]*_down_bound_U[_nbr_elements_1D-1][1] + _down_bound_U[_nbr_elements_1D-1][2]*_down_bound_U[_nbr_elements_1D-1][2])/_down_bound_U[_nbr_elements_1D-1][0];
  sound_speed = sqrt(1.4*_down_bound_Pi[_nbr_elements_1D-1]/_down_bound_U[_nbr_elements_1D-1][0]);
  file << " " << velocity_mag;
  file << " " << velocity_mag/sound_speed;
  file << std::endl;



  for (int line = 0; line < _nbr_elements_1D; line++)
  {
    file << -0.5*_delta_s << " " << (line + 0.5)*_delta_s;
    file << " " << _left_bound_U[line][0];
    file << " " << _left_bound_Pi[line];
    velocity_mag = sqrt(_left_bound_U[line][1]*_left_bound_U[line][1] + _left_bound_U[line][2]*_left_bound_U[line][2])/_left_bound_U[line][0];
    sound_speed = sqrt(1.4*_left_bound_Pi[line]/_left_bound_U[line][0]);
    file << " " << velocity_mag;
    file << " " << velocity_mag/sound_speed;
    file << std::endl;

    for (int column = 0; column < _nbr_elements_1D; column++)
    {
      file << (column + 0.5)*_delta_s << " " << (line + 0.5)*_delta_s;
      file << " " << _U[line][column][0];
      file << " " << _Pi[line][column];
      velocity_mag = sqrt(_U[line][column][1]*_U[line][column][1] + _U[line][column][2]*_U[line][column][2])/_U[line][column][0];
      sound_speed = sqrt(1.4*_Pi[line][column]/_U[line][column][0]);
      file << " " << velocity_mag;
      file << " " << velocity_mag/sound_speed;
      file << " " << _a_LR[line][column] << " " << _u_star[line][column] << " " << _Pi_star_LR[line][column] << " " << _L[line][column];
      file << std::endl;
    }

    file << 1.0 + 0.5*_delta_s << " " << (line + 0.5)*_delta_s;
    file << " " << _right_bound_U[line][0];
    file << " " << _right_bound_Pi[line];
    velocity_mag = sqrt(_right_bound_U[line][1]*_right_bound_U[line][1] + _right_bound_U[line][2]*_right_bound_U[line][2])/_right_bound_U[line][0];
    sound_speed = sqrt(1.4*_right_bound_Pi[line]/_right_bound_U[line][0]);
    file << " " << velocity_mag;
    file << " " << velocity_mag/sound_speed;
    file << std::endl;
  }



  file << -0.5*_delta_s << " " << 1.0 + 0.5*_delta_s;
  file << " " << _up_bound_U[0][0]; //Le coin a la même valeur que le premier élément du haut
  file << " " << _up_bound_Pi[0];
  velocity_mag = sqrt(_up_bound_U[0][1]*_up_bound_U[0][1] + _up_bound_U[0][2]*_up_bound_U[0][2])/_up_bound_U[0][0];
  sound_speed = sqrt(1.4*_up_bound_Pi[0]/_up_bound_U[0][0]);
  file << " " << velocity_mag;
  file << " " << velocity_mag/sound_speed;
  file << std::endl;

  for (int column = 0; column < _nbr_elements_1D; column++)
  {
    file << (column + 0.5)*_delta_s << " " << 1.0 + 0.5*_delta_s;
    file << " " << _up_bound_U[column][0];
    file << " " << _up_bound_Pi[column];
    velocity_mag = sqrt(_up_bound_U[column][1]*_up_bound_U[column][1] + _up_bound_U[column][2]*_up_bound_U[column][2])/_up_bound_U[column][0];
    sound_speed = sqrt(1.4*_up_bound_Pi[column]/_up_bound_U[column][0]);
    file << " " << velocity_mag;
    file << " " << velocity_mag/sound_speed;
    file << std::endl;
  }

  file << 1.0 + 0.5*_delta_s << " " << 1.0 + 0.5*_delta_s;
  file << " " << _up_bound_U[_nbr_elements_1D-1][0]; //Le coin a la même valeur que le dernier élément du haut
  file << " " << _up_bound_Pi[_nbr_elements_1D-1];
  velocity_mag = sqrt(_up_bound_U[_nbr_elements_1D-1][1]*_up_bound_U[_nbr_elements_1D-1][1] + _up_bound_U[_nbr_elements_1D-1][2]*_up_bound_U[_nbr_elements_1D-1][2])/_up_bound_U[_nbr_elements_1D-1][0];
  sound_speed = sqrt(1.4*_up_bound_Pi[_nbr_elements_1D-1]/_up_bound_U[_nbr_elements_1D-1][0]);
  file << " " << velocity_mag;
  file << " " << velocity_mag/sound_speed;
  file << std::endl;


  // file << "# état du système à t = " << _time << std::endl;
  // file << -0.5*_delta_s << " " << -0.5*_delta_s;
  // for (int iVar = 0; iVar < 4; iVar++)
  // {
  //   file << " " << _down_bound_U[0][iVar]; //Le coin a la même valeur que le premier élément du bas
  // }
  // file << " " << _down_bound_Pi[0] << std::endl;
  //
  // for (int column = 0; column < _nbr_elements_1D; column++)
  // {
  //   file << (column + 0.5)*_delta_s << " " << -0.5*_delta_s;
  //   for (int iVar = 0; iVar < 4; iVar++)
  //   {
  //     file << " " << _down_bound_U[column][iVar];
  //   }
  //   file << " " << _down_bound_Pi[column] << std::endl;
  // }
  //
  // file << 1.0 + 0.5*_delta_s << " " << -0.5*_delta_s;
  // for (int iVar = 0; iVar < 4; iVar++)
  // {
  //   file << " " << _down_bound_U[_nbr_elements_1D-1][iVar]; //Le coin a la même valeur que le dernier élément du bas
  // }
  // file << " " << _down_bound_Pi[_nbr_elements_1D-1] << std::endl;
  //
  //
  //
  // for (int line = 0; line < _nbr_elements_1D; line++)
  // {
  //   file << -0.5*_delta_s << " " << (line + 0.5)*_delta_s;
  //   for (int iVar = 0; iVar < 4; iVar++)
  //   {
  //     file << " " << _left_bound_U[line][iVar];
  //   }
  //   file << " " << _left_bound_Pi[line] << std::endl;
  //
  //   for (int column = 0; column < _nbr_elements_1D; column++)
  //   {
  //     file << (column + 0.5)*_delta_s << " " << (line + 0.5)*_delta_s;
  //     for (int iVar = 0; iVar < 4; iVar++)
  //     {
  //       file << " " << _U[line][column][iVar];
  //     }
  //     file << " " << _Pi[line][column] << std::endl;
  //   }
  //
  //   file << 1.0 + 0.5*_delta_s << " " << (line + 0.5)*_delta_s;
  //   for (int iVar = 0; iVar < 4; iVar++)
  //   {
  //     file << " " << _right_bound_U[line][iVar];
  //   }
  //   file << " " << _right_bound_Pi[line] << std::endl;
  // }
  //
  //
  //
  // file << -0.5*_delta_s << " " << 1.0 + 0.5*_delta_s;
  // for (int iVar = 0; iVar < 4; iVar++)
  // {
  //   file << " " << _up_bound_U[0][iVar]; //Le coin a la même valeur que le premier élément du haut
  // }
  // file << " " << _up_bound_Pi[0] << std::endl;
  //
  // for (int column = 0; column < _nbr_elements_1D; column++)
  // {
  //   file << (column + 0.5)*_delta_s << " " << 1.0 + 0.5*_delta_s;
  //   for (int iVar = 0; iVar < 4; iVar++)
  //   {
  //     file << " " << _up_bound_U[column][iVar];
  //   }
  //   file << " " << _up_bound_Pi[column] << std::endl;
  // }
  //
  // file << 1.0 + 0.5*_delta_s << " " << 1.0 + 0.5*_delta_s;
  // for (int iVar = 0; iVar < 4; iVar++)
  // {
  //   file << " " << _up_bound_U[_nbr_elements_1D-1][iVar]; //Le coin a la même valeur que le dernier élément du haut
  // }
  // file << " " << _up_bound_Pi[_nbr_elements_1D-1] << std::endl;

  file.close();
}

void Probleme2D::TimeIteration(int time_it)
{
  Update_a_u_v_star();
  Update_theta();
  Update_Pi_star();

  //Evaluation de la cfl pour l'acoustic step
  double max_for_cfl = 0.0;
  double buffer;
  for (int line = 0; line < _nbr_elements_1D; line++)
  {
    for (int column = 0; column < _nbr_elements_1D; column++)
    {
      buffer = std::max(std::max(_a_LR[line][column],_a_LR[line][column+1]),std::max(_a_DU[line][column],_a_DU[line+1][column]))/_U[line][column][0];
      if (buffer > max_for_cfl)
      {
        max_for_cfl = buffer;
      }
    }
  }

  _Dt_on_Ds = _cfl/(2*max_for_cfl);

  //Evaluation de la cfl pour la transport step
  max_for_cfl = 0.0;

  for (int line = 0; line < _nbr_elements_1D; line++)
  {
    for (int column = 0; column < _nbr_elements_1D; column++)
    {
      buffer = abs(_u_star[line][column]) + abs(_u_star[line][column+1]) + abs(_v_star[line][column]) + abs(_v_star[line+1][column]);
      if (buffer > max_for_cfl)
      {
        max_for_cfl = buffer;
      }
    }
  }

  _Dt_on_Ds = std::min(_cfl/max_for_cfl,_Dt_on_Ds);

  _time += _Dt_on_Ds*_delta_s;

  Update_L();

  AcousticStep();

  Update_Phi_interface();

  TransportStep();

  Update_u_v_Pi();

  Update_CL();

  SaveIteration(time_it);

}

void Probleme2D::Solve()
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
double PressureToRhoE2D(double p, double rho, double rho_u, double rho_v)
{
  double rho_E = 2.5*p + 0.5*(rho_u*rho_u + rho_v*rho_v)/rho;
  return rho_E;
}
