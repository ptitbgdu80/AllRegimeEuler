#include "AREuler.h"

//Constructeur
Problem2D::Problem2D(int nbr_elements_1D, double t_final, int choice_theta, std::string file_name, int choice_test_case, int choice_solver)
{
  _nbr_elements_1D = nbr_elements_1D;
  _t_final = t_final;
  _choice_theta = choice_theta;
  _file_name = file_name;
  double pi = 3.141592653589793;
  if (choice_test_case == 2)
  {
    _length = 2*pi;
  }
  else
  {
    _length = 1.0;
  }
  _delta_s = _length/nbr_elements_1D;
  _choice_test_case = choice_test_case;
  _choice_solver = choice_solver;

  _cfl = 0.95;

  _u.resize(_nbr_elements_1D);
  _v.resize(_nbr_elements_1D);
  _Pi.resize(_nbr_elements_1D); // Pi contient les pressions calculées à partir de U
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

  switch (choice_solver)
  {
    case 0 :

    _L.resize(nbr_elements_1D); // L[j] = 1 + Dt/Ds*somme(u*jk);

    _a_LR.resize(nbr_elements_1D);
    _theta_LR.resize(nbr_elements_1D);
    _u_star.resize(nbr_elements_1D);
    _Pi_star_LR.resize(nbr_elements_1D);
    _Phi_interface_LR.resize(nbr_elements_1D);
    for (int line = 0; line < nbr_elements_1D; line++)
    {
      _L[line].resize(nbr_elements_1D);

      _a_LR[line].resize(nbr_elements_1D+1);
      _theta_LR[line].resize(nbr_elements_1D+1);
      _u_star[line].resize(nbr_elements_1D+1);
      _Pi_star_LR[line].resize(nbr_elements_1D+1);
      _Phi_interface_LR[line].resize(nbr_elements_1D+1);
      for (int column = 0; column < nbr_elements_1D+1; column++)
      {
        _Phi_interface_LR[line][column].resize(4);
      }
    }

    _a_DU.resize(nbr_elements_1D+1);
    _theta_DU.resize(nbr_elements_1D+1);
    _v_star.resize(nbr_elements_1D+1);
    _Pi_star_DU.resize(nbr_elements_1D+1);
    _Phi_interface_DU.resize(nbr_elements_1D+1);
    for (int line = 0; line < nbr_elements_1D+1; line++)
    {
      _a_DU[line].resize(nbr_elements_1D);
      _theta_DU[line].resize(nbr_elements_1D);
      _v_star[line].resize(nbr_elements_1D);
      _Pi_star_DU[line].resize(nbr_elements_1D);
      _Phi_interface_DU[line].resize(nbr_elements_1D);
      for (int column = 0; column < nbr_elements_1D; column++)
      {
        _Phi_interface_DU[line][column].resize(4);
      }
    }
    break;

    case 1 :

    _Fx.resize(nbr_elements_1D);
    _Fy.resize(nbr_elements_1D);
    _left_bound_Fx.resize(nbr_elements_1D);
    _right_bound_Fx.resize(nbr_elements_1D);
    _down_bound_Fy.resize(nbr_elements_1D);
    _up_bound_Fy.resize(nbr_elements_1D);

    for (int line = 0; line < nbr_elements_1D; line++)
    {
      _Fx[line].resize(nbr_elements_1D);
      _Fy[line].resize(nbr_elements_1D);

      _left_bound_Fx[line].resize(4);
      _right_bound_Fx[line].resize(4);
      _down_bound_Fy[line].resize(4); // Devrait être mis dans une boucle sur les colonnes...
      _up_bound_Fy[line].resize(4); // idem

      for (int column = 0; column < nbr_elements_1D; column++)
      {
        _Fx[line][column].resize(4);
        _Fy[line][column].resize(4);
      }
    }

    _max_c_LR.resize(nbr_elements_1D);
    _flow_LR.resize(nbr_elements_1D);

    for (int line = 0; line < nbr_elements_1D; line++)
    {
      _max_c_LR[line].resize(nbr_elements_1D+1);
      _flow_LR[line].resize(nbr_elements_1D+1);
      for (int column = 0; column < nbr_elements_1D+1; column++)
      {
        _flow_LR[line][column].resize(4);
      }
    }

    _max_c_DU.resize(nbr_elements_1D+1);
    _flow_DU.resize(nbr_elements_1D+1);

    for (int line = 0; line < nbr_elements_1D+1; line++)
    {
      _max_c_DU[line].resize(nbr_elements_1D);
      _flow_DU[line].resize(nbr_elements_1D);
      for (int column = 0; column < nbr_elements_1D; column++)
      {
        _flow_DU[line][column].resize(4);
      }
    }
    break;

    default :
    std::cout << "Le choix du solveur doit être 0 ou 1" << std::endl;
    exit(1);
  }



  //Conditions initiales
  switch (choice_test_case)
  {
    case 0 :
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
        _u[line][column] = u;
        _v[line][column] = v;
        _Pi[line][column] = p;

        _U[line][column][0] = rho;
        _U[line][column][1] = rho*u;
        _U[line][column][2] = rho*v;
        _U[line][column][3] = PressureToRhoE2D(p, _U[line][column][0], _U[line][column][1], _U[line][column][2]);
      }
    }
    break;

    case 1 :
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
    break;

    case 2 :
    for (int line = 0; line < nbr_elements_1D; line++)
    {
      for (int column = 0; column < nbr_elements_1D; column++)
      {
        double x ,y ,rho, u, v, p;
        x = column*_delta_s;
        y = line*_delta_s;

        rho = 1.96;
        u = -sin(y);
        v = sin(x);
        p = 1.4;

        _u[line][column] = u;
        _v[line][column] = v;
        _Pi[line][column] = p;

        _U[line][column][0] = rho;
        _U[line][column][1] = rho*u;
        _U[line][column][2] = rho*v;
        _U[line][column][3] = PressureToRhoE2D(p, _U[line][column][0], _U[line][column][1], _U[line][column][2]);
      }
    }
    break;

    case 3 :
    for (int line = 0; line < nbr_elements_1D; line++)
    {
      for (int column = 0; column < nbr_elements_1D; column++)
      {
        double x ,rho, u, v, p;
        x = (column+0.5)*_delta_s;

        if (x < 0.5)
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

        _u[line][column] = u;
        _v[line][column] = v;
        _Pi[line][column] = p;

        _U[line][column][0] = rho;
        _U[line][column][1] = rho*u;
        _U[line][column][2] = rho*v;
        _U[line][column][3] = PressureToRhoE2D(p, _U[line][column][0], _U[line][column][1], _U[line][column][2]);
      }
    }
    break;


    default :
    std::cout << "Le choix du cas test doit être 0, 1, ou 2" << std::endl;
    exit(1);
  }

  //Conditions aux bords
  Update_CL();

  //Initialisation des theta
  switch (choice_theta)
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
Problem2D::~Problem2D()
{}



//Methodes de la classe
void Problem2D::Update_u_v_Pi()
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

void Problem2D::Update_a_u_v_star()
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

void Problem2D::Update_theta()
{
  double cL, cR, cD, cU;
  switch (_choice_theta)
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

void Problem2D::Update_Pi_star()
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

void Problem2D::Update_L()
{
  for (int line = 0; line < _nbr_elements_1D; line++)
  {
    for (int column = 0; column < _nbr_elements_1D; column++)
    {
      _L[line][column] = 1.0 + _Dt_on_Ds*(_u_star[line][column+1] - _u_star[line][column] + _v_star[line+1][column] - _v_star[line][column]);
    }
  }
}

void Problem2D::Update_Phi_interface()
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

void Problem2D::Update_CL()
{
  switch (_choice_test_case)
  {
    case 0 :
    for (int index = 0; index < _nbr_elements_1D; index++)
    {
      _left_bound_Pi[index] = _Pi[index][0];
      _right_bound_Pi[index] = _Pi[index][_nbr_elements_1D-1];
      _down_bound_Pi[index] = _Pi[0][index];
      _up_bound_Pi[index] = _Pi[_nbr_elements_1D-1][index];

      _left_bound_U[index][0] = _U[index][0][0];
      _left_bound_U[index][1] = 0.0;
      _left_bound_U[index][2] = 0.0;
      _left_bound_U[index][3] = PressureToRhoE(_left_bound_Pi[index], _left_bound_U[index][0], 0, 0);

      _right_bound_U[index][0] = _U[index][_nbr_elements_1D-1][0];
      _right_bound_U[index][1] = 0.0;
      _right_bound_U[index][2] = 0.0;
      _right_bound_U[index][3] = PressureToRhoE(_right_bound_Pi[index], _right_bound_U[index][0], 0, 0);

      _down_bound_U[index][0] = _U[0][index][0];
      _down_bound_U[index][1] = 0.0;
      _down_bound_U[index][2] = 0.0;
      _down_bound_U[index][3] = PressureToRhoE(_down_bound_Pi[index], _down_bound_U[index][0], 0, 0);

      _up_bound_U[index][0] = _U[_nbr_elements_1D-1][index][0];
      _up_bound_U[index][1] = 0.0;
      _up_bound_U[index][2] = 0.0;
      _up_bound_U[index][3] = PressureToRhoE(_up_bound_Pi[index], _up_bound_U[index][0], 0, 0);

      _left_bound_u[index] = 0.0;
      _right_bound_u[index] = 0.0;
      _down_bound_v[index] = 0.0;
      _up_bound_v[index] = 0.0;
    }
    break;

    case 1 :
    case 3 :
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
    break;

    case 2 :
    for (int index = 0; index < _nbr_elements_1D; index++)
    {
      // CL périodiques
      _left_bound_Pi[index] = _Pi[index][_nbr_elements_1D-1];
      _right_bound_Pi[index] = _Pi[index][0];
      _down_bound_Pi[index] = _Pi[_nbr_elements_1D-1][index];
      _up_bound_Pi[index] = _Pi[0][index];

      _left_bound_U[index][0] = _U[index][_nbr_elements_1D-1][0];
      _left_bound_U[index][1] = _U[index][_nbr_elements_1D-1][1];
      _left_bound_U[index][2] = _U[index][_nbr_elements_1D-1][2];
      _left_bound_U[index][3] = PressureToRhoE(_left_bound_Pi[index], _left_bound_U[index][0], _left_bound_U[index][1], _left_bound_U[index][2]);

      _right_bound_U[index][0] = _U[index][0][0];
      _right_bound_U[index][1] = _U[index][0][1];
      _right_bound_U[index][2] = _U[index][0][2];
      _right_bound_U[index][3] = PressureToRhoE(_right_bound_Pi[index], _right_bound_U[index][0], _right_bound_U[index][1], _right_bound_U[index][2]);

      _down_bound_U[index][0] = _U[_nbr_elements_1D-1][index][0];
      _down_bound_U[index][1] = _U[_nbr_elements_1D-1][index][1];
      _down_bound_U[index][2] = _U[_nbr_elements_1D-1][index][2];
      _down_bound_U[index][3] = PressureToRhoE(_down_bound_Pi[index], _down_bound_U[index][0], _down_bound_U[index][1], _down_bound_U[index][2]);

      _up_bound_U[index][0] = _U[0][index][0];
      _up_bound_U[index][1] = _U[0][index][1];
      _up_bound_U[index][2] = _U[0][index][2];
      _up_bound_U[index][3] = PressureToRhoE(_up_bound_Pi[index], _up_bound_U[index][0], _up_bound_U[index][1], _up_bound_U[index][2]);

      _left_bound_u[index] = _u[index][_nbr_elements_1D-1];
      _right_bound_u[index] = _u[index][0];
      _down_bound_v[index] = _v[_nbr_elements_1D-1][index];
      _up_bound_v[index] = _v[0][index];
    }
    break;

    default :
    std::cout << "Le choix du cas test doit être 0, 1, ou 2" << std::endl;
    exit(1);
  }

  if (_choice_solver == 1)
  {
    for (int index = 0; index < _nbr_elements_1D; index++)
    {
      _left_bound_Fx[index][0] = _left_bound_U[index][1];
      _left_bound_Fx[index][1] = _left_bound_U[index][1]*_left_bound_U[index][1]/_left_bound_U[index][0] + _left_bound_Pi[index];
      _left_bound_Fx[index][2] = _left_bound_U[index][1]*_left_bound_U[index][2]/_left_bound_U[index][0];
      _left_bound_Fx[index][3] = _left_bound_U[index][1]/_left_bound_U[index][0]*(_left_bound_U[index][3] + _left_bound_Pi[index]);

      _right_bound_Fx[index][0] = _right_bound_U[index][1];
      _right_bound_Fx[index][1] = _right_bound_U[index][1]*_right_bound_U[index][1]/_right_bound_U[index][0] + _right_bound_Pi[index];
      _right_bound_Fx[index][2] = _right_bound_U[index][1]*_right_bound_U[index][2]/_right_bound_U[index][0];
      _right_bound_Fx[index][3] = _right_bound_U[index][1]/_right_bound_U[index][0]*(_right_bound_U[index][3] + _right_bound_Pi[index]);

      _down_bound_Fy[index][0] = _down_bound_U[index][2];
      _down_bound_Fy[index][1] = _down_bound_U[index][1]*_down_bound_U[index][2]/_down_bound_U[index][0];
      _down_bound_Fy[index][2] = _down_bound_U[index][2]*_down_bound_U[index][2]/_down_bound_U[index][0] + _down_bound_Pi[index];
      _down_bound_Fy[index][3] = _down_bound_U[index][2]/_down_bound_U[index][0]*(_down_bound_U[index][3] + _down_bound_Pi[index]);

      _up_bound_Fy[index][0] = _up_bound_U[index][2];
      _up_bound_Fy[index][1] = _up_bound_U[index][1]*_up_bound_U[index][2]/_up_bound_U[index][0];
      _up_bound_Fy[index][2] = _up_bound_U[index][2]*_up_bound_U[index][2]/_up_bound_U[index][0] + _up_bound_Pi[index];
      _up_bound_Fy[index][3] = _up_bound_U[index][2]/_up_bound_U[index][0]*(_up_bound_U[index][3] + _up_bound_Pi[index]);
    }
  }
}

void Problem2D::AcousticStep()
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

void Problem2D::TransportStep()
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

void Problem2D::Update_F()
{
  for (int line = 0; line < _nbr_elements_1D; line++)
  {
    for (int column = 0; column < _nbr_elements_1D; column++)
    {
      _Fx[line][column][0] = _U[line][column][1];
      _Fx[line][column][1] = _U[line][column][1]*_U[line][column][1]/_U[line][column][0] + _Pi[line][column];
      _Fx[line][column][2] = _U[line][column][1]*_U[line][column][2]/_U[line][column][0];
      _Fx[line][column][3] = _U[line][column][1]/_U[line][column][0]*(_U[line][column][3] + _Pi[line][column]);

      _Fy[line][column][0] = _U[line][column][2];
      _Fy[line][column][1] = _U[line][column][1]*_U[line][column][2]/_U[line][column][0];
      _Fy[line][column][2] = _U[line][column][2]*_U[line][column][2]/_U[line][column][0] + _Pi[line][column];
      _Fy[line][column][3] = _U[line][column][2]/_U[line][column][0]*(_U[line][column][3] + _Pi[line][column]);
    }
  }
}

void Problem2D::Update_max_c()
{
  for (int line = 0; line < _nbr_elements_1D; line++)
  {
    double cL = sqrt(1.4*_left_bound_Pi[line]/_left_bound_U[line][0]);
    double cR = sqrt(1.4*_Pi[line][0]/_U[line][0][0]);
    _max_c_LR[line][0] = std::max(cL,cR);

    for (int column = 1; column < _nbr_elements_1D; column++)
    {
      cL = cR;
      cR = sqrt(1.4*_Pi[line][column]/_U[line][column][0]);
      _max_c_LR[line][column] = std::max(cL, cR);
    }

    cL = cR;
    cR = sqrt(1.4*_right_bound_Pi[line]/_right_bound_U[line][0]);
    _max_c_LR[line][_nbr_elements_1D] = std::max(cL, cR);
  }


  for (int column = 0; column < _nbr_elements_1D; column++)
  {
    double cD = sqrt(1.4*_down_bound_Pi[column]/_down_bound_U[column][0]);
    double cU = sqrt(1.4*_Pi[0][column]/_U[0][column][0]);
    _max_c_DU[0][column] = std::max(cD, cU);

    for (int line = 1; line < _nbr_elements_1D; line++)
    {
      cD = cU;
      cU = sqrt(1.4*_Pi[line][column]/_U[line][column][0]);
      _max_c_DU[line][column] = std::max(cD, cU);
    }

    cD = cU;
    cU = sqrt(1.4*_up_bound_Pi[column]/_up_bound_U[column][0]);
    _max_c_DU[_nbr_elements_1D][column] = std::max(cD, cU);
  }
}

void Problem2D::Update_flow()
{
  for (int iVar = 0; iVar < 4; iVar++)
  {
    for (int line = 0; line < _nbr_elements_1D; line++)
    {
      _flow_LR[line][0][iVar] = 0.5*(_Fx[line][0][iVar] + _left_bound_Fx[line][iVar] - _max_c_LR[line][iVar]*(_U[line][0][iVar] - _left_bound_U[line][iVar]));
      for (int column = 1; column < _nbr_elements_1D; column++)
      {
        _flow_LR[line][column][iVar] = 0.5*(_Fx[line][column][iVar] + _Fx[line][column-1][iVar] - _max_c_LR[line][column]*(_U[line][column][iVar] - _U[line][column-1][iVar]));
      }
      _flow_LR[line][_nbr_elements_1D][iVar] = 0.5*(_right_bound_Fx[line][iVar] + _Fx[line][_nbr_elements_1D-1][iVar] - _max_c_LR[line][_nbr_elements_1D]*(_right_bound_U[line][iVar] - _U[line][_nbr_elements_1D-1][iVar]));
    }

    for (int column = 0; column < _nbr_elements_1D; column++)
    {
      _flow_DU[0][column][iVar] = 0.5*(_Fy[0][column][iVar] + _down_bound_Fy[column][iVar] - _max_c_DU[0][column]*(_U[0][column][iVar] - _down_bound_U[column][iVar]));
      for (int line = 1; line < _nbr_elements_1D; line++)
      {
        _flow_DU[line][column][iVar] = 0.5*(_Fy[line][column][iVar] + _Fy[line-1][column][iVar] - _max_c_DU[line][column]*(_U[line][column][iVar] - _U[line-1][column][iVar]));
      }
      _flow_DU[_nbr_elements_1D][column][iVar] = 0.5*(_up_bound_Fy[column][iVar] + _Fy[_nbr_elements_1D-1][column][iVar] - _max_c_DU[_nbr_elements_1D][column]*(_up_bound_U[column][iVar] - _U[_nbr_elements_1D-1][column][iVar]));
    }
  }
}

void Problem2D::RusanovStep()
{
  for (int iVar = 0; iVar < 4; iVar++)
  {
    for (int line = 0; line < _nbr_elements_1D; line++)
    {
      for (int column = 0; column < _nbr_elements_1D; column++)
      {
        _U[line][column][iVar] -= _Dt_on_Ds*(_flow_LR[line][column+1][iVar] - _flow_LR[line][column][iVar] + _flow_DU[line+1][column][iVar] - _flow_DU[line][column][iVar]);
      }
    }
  }
}

void Problem2D::SaveIteration(int time_it)
{
  if (time_it == 0)
  {
    system(("rm -rf " + _file_name).c_str());
    system(("mkdir " + _file_name).c_str());
  }

  std::ofstream file;
  file.open(_file_name + "/it_" + std::to_string(time_it), std::ios::out);

  file << "# état du système à t = " << _time << std::endl;
  if (_choice_test_case == 3)
  {
    int line = _nbr_elements_1D/2;
    file << -0.5*_delta_s;
    file << " " << _left_bound_U[line][0];
    file << " " << _left_bound_Pi[line];
    double velocity_mag = sqrt(_left_bound_U[line][1]*_left_bound_U[line][1] + _left_bound_U[line][2]*_left_bound_U[line][2])/_left_bound_U[line][0];
    double sound_speed = sqrt(1.4*_left_bound_Pi[line]/_left_bound_U[line][0]);
    file << " " << velocity_mag;
    file << " " << velocity_mag/sound_speed;
    file << std::endl;

    for (int column = 0; column < _nbr_elements_1D; column++)
    {
      file << (column + 0.5)*_delta_s;
      file << " " << _U[line][column][0];
      file << " " << _Pi[line][column];
      velocity_mag = sqrt(_U[line][column][1]*_U[line][column][1] + _U[line][column][2]*_U[line][column][2])/_U[line][column][0];
      sound_speed = sqrt(1.4*_Pi[line][column]/_U[line][column][0]);
      file << " " << velocity_mag;
      file << " " << velocity_mag/sound_speed;
      file << std::endl;
    }

    file << _length + 0.5*_delta_s;
    file << " " << _right_bound_U[line][0];
    file << " " << _right_bound_Pi[line];
    velocity_mag = sqrt(_right_bound_U[line][1]*_right_bound_U[line][1] + _right_bound_U[line][2]*_right_bound_U[line][2])/_right_bound_U[line][0];
    sound_speed = sqrt(1.4*_right_bound_Pi[line]/_right_bound_U[line][0]);
    file << " " << velocity_mag;
    file << " " << velocity_mag/sound_speed;
    file << std::endl;
  }
  else
  {
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

    file << _length + 0.5*_delta_s << " " << -0.5*_delta_s;
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
        file << std::endl;
      }

      file << _length + 0.5*_delta_s << " " << (line + 0.5)*_delta_s;
      file << " " << _right_bound_U[line][0];
      file << " " << _right_bound_Pi[line];
      velocity_mag = sqrt(_right_bound_U[line][1]*_right_bound_U[line][1] + _right_bound_U[line][2]*_right_bound_U[line][2])/_right_bound_U[line][0];
      sound_speed = sqrt(1.4*_right_bound_Pi[line]/_right_bound_U[line][0]);
      file << " " << velocity_mag;
      file << " " << velocity_mag/sound_speed;
      file << std::endl;
    }



    file << -0.5*_delta_s << " " << _length + 0.5*_delta_s;
    file << " " << _up_bound_U[0][0]; //Le coin a la même valeur que le premier élément du haut
    file << " " << _up_bound_Pi[0];
    velocity_mag = sqrt(_up_bound_U[0][1]*_up_bound_U[0][1] + _up_bound_U[0][2]*_up_bound_U[0][2])/_up_bound_U[0][0];
    sound_speed = sqrt(1.4*_up_bound_Pi[0]/_up_bound_U[0][0]);
    file << " " << velocity_mag;
    file << " " << velocity_mag/sound_speed;
    file << std::endl;

    for (int column = 0; column < _nbr_elements_1D; column++)
    {
      file << (column + 0.5)*_delta_s << " " << _length + 0.5*_delta_s;
      file << " " << _up_bound_U[column][0];
      file << " " << _up_bound_Pi[column];
      velocity_mag = sqrt(_up_bound_U[column][1]*_up_bound_U[column][1] + _up_bound_U[column][2]*_up_bound_U[column][2])/_up_bound_U[column][0];
      sound_speed = sqrt(1.4*_up_bound_Pi[column]/_up_bound_U[column][0]);
      file << " " << velocity_mag;
      file << " " << velocity_mag/sound_speed;
      file << std::endl;
    }

    file << _length + 0.5*_delta_s << " " << _length + 0.5*_delta_s;
    file << " " << _up_bound_U[_nbr_elements_1D-1][0]; //Le coin a la même valeur que le dernier élément du haut
    file << " " << _up_bound_Pi[_nbr_elements_1D-1];
    velocity_mag = sqrt(_up_bound_U[_nbr_elements_1D-1][1]*_up_bound_U[_nbr_elements_1D-1][1] + _up_bound_U[_nbr_elements_1D-1][2]*_up_bound_U[_nbr_elements_1D-1][2])/_up_bound_U[_nbr_elements_1D-1][0];
    sound_speed = sqrt(1.4*_up_bound_Pi[_nbr_elements_1D-1]/_up_bound_U[_nbr_elements_1D-1][0]);
    file << " " << velocity_mag;
    file << " " << velocity_mag/sound_speed;
    file << std::endl;
  }

  file.close();
}

void Problem2D::TimeIteration(int time_it)
{
  double max_for_cfl, buffer;
  switch (_choice_solver)
  {
    case 0 :
    Update_a_u_v_star();
    Update_theta();
    Update_Pi_star();

    //Evaluation de la cfl pour l'acoustic step
    max_for_cfl = 0.0;
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

    _Dt_on_Ds = _cfl/(8*max_for_cfl);

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

    _Dt_on_Ds = std::min(_cfl/(4*max_for_cfl),_Dt_on_Ds);
    _time += _Dt_on_Ds*_delta_s;

    Update_L();
    AcousticStep();
    Update_Phi_interface();
    TransportStep();
    break;

    case 1 :
    Update_F();
    Update_max_c();
    Update_flow();

    //Evaluation de la cfl
    max_for_cfl = 0.0;
    for (int line = 0; line < _nbr_elements_1D; line++)
    {
      for (int column = 0; column < _nbr_elements_1D; column++)
      {
        buffer = std::max(_max_c_LR[line][column],_max_c_LR[line][column+1]),std::max(_max_c_DU[line][column],_max_c_DU[line+1][column]);
        if (buffer > max_for_cfl)
        {
          max_for_cfl = buffer;
        }
      }
    }

    _Dt_on_Ds = _cfl/(8*max_for_cfl);
    _time += _Dt_on_Ds*_delta_s;

    RusanovStep();
    break;

    default :
    std::cout << "le choix du solveur doit être 0 ou 1" << std::endl;
    exit(1);
  }

  Update_u_v_Pi();
  Update_CL();
  if (time_it%20 == 0)
  {
    SaveIteration(time_it);
  }
}

void Problem2D::Solve()
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
