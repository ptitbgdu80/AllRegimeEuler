#include "AREuler.h"

//Constructeurs
Probleme1D::Probleme1D()
{}

Probleme1D::Probleme1D(int nbr_elements, double delta_t, double t_final, int flow_choice, std::string file_name)
{
  //Cette classe contient 3 variables conservatives rho, rho*u et rho*E

  _nbr_elements = nbr_elements;
  _delta_t = delta_t;
  _t_final = t_final;
  _flow_choice = flow_choice;
  _file_name = file_name;
  _delta_x = 1./nbr_elements;

  _U.resize(nbr_elements);
  _W.resize(nbr_elements);
  _flows.resize(nbr_elements);
  for (int elem_i = 0; elem_i < nbr_elements; elem_i++)
  {
    _U[elem_i].resize(4);
    _W[elem_i].resize(5);
    for (int jVar = 0; jVar < 5; jVar++)
    {
      _U[elem_i][jVar] = 0;
    }
  }

}



//Destructeur
Probleme1D::~Probleme1D()
{}



//Methodes générales
void Probleme1D::AcousticStep()
{


  //Mise à jour de U à partir de W
  for (int elem_j = 0; elem_j < _nbr_elements; elem_j++)
  {
    double inverse_Lj = 1./(1. - coeff*flow_contribution[elem_j][0]);

  }

}

void Probleme1D::SaveIteration(int time_it)
{

}

std::vector<double> Probleme1D::LeftBoundAcousticFlow()
{
  double rho, u, E;
  rho = 1.2; //masse volumique de l'air
  u = 10.;
  E = 1.;
  return {rho, rho*u, rho*E};
}

std::vector<double> Probleme1D::RightBoundAcousticFlow()
{
  double rho, u, E;
  rho = 1.2; //masse volumique de l'air
  u = 10.;
  E = 1.;
  return {rho, rho*u, rho*E};
}

void Probleme1D::ClassicFVMainLoop()
{
  int nbr_iterations = floor(_t_final/_delta_t + 1);
  for (int time_it = 0; time_it < nbr_iterations; time_it++)
  {
    Probleme1D::UpdateFlows();
    Probleme1D::AcousticStep();
  }
}
