#include "AREuler.h"

//Constructeurs
Probleme1D::Probleme1D()
{}

TroisVariables::TroisVariables(int nbr_elements, double delta_t, int choix_flux)
{
  //Cette classe contient 3 variables conservatives rho, rho*u et rho*E

  _nbr_elements = nbr_elements;
  _delta_t = delta_t;
  _choix_flux = choix_flux;
  _delta_x = 1./nbr_elements;
  _nbr_variables = 3;

  _U.resize(nbr_elements);
  _flows.resize(nbr_elements);
  for (int elem_i = 0; elem_i < nbr_elements; elem_i++)
  {
    _U[elem_i].resize(_nbr_variables);
    _flows[elem_i].resize(_nbr_variables);
    for (int jVar = 0; jVar < _nbr_variables; jVar++)
    _U[elem_i][jVar] = 0;
  }

}



//Destructeur
Probleme1D::~Probleme1D()
{}



//Methodes générales
void Probleme1D::ClassicGodunovIteration()
{
  //Initialisation d'un vecteur qui contient les quantités (F_i+1/2 - F_i-1/2)
  std::vector<std::vector<double>> flow_contribution;
  flow_contribution.resize(_nbr_elements);

  double _nbr_variables = _U[0].size();
  for (int elem_i = 0; elem_i < _nbr_elements; elem_i++)
  {
    flow_contribution[elem_i].resize(_nbr_variables);
    for (int jVar = 0; jVar < _nbr_variables; jVar++)
    {
      flow_contribution[elem_i][jVar] = 0;
    }
  }

  //Somme des contributions avec un parcours par interfaces
  for (int elem_i = 1; elem_i < _nbr_elements; elem_i++)
  {
    std::vector<double> flow = RightInterfaceFlow(elem_i);
    for (int jVar = 0; jVar < _nbr_variables; jVar++)
    {
      flow_contribution[elem_i-1][jVar] -= flow[jVar];
      flow_contribution[elem_i][jVar] += flow[jVar];
    }
  }

  //Contribution des bords
  std::vector<double> left_flow, right_flow;
  left_flow = LeftBoundFlow();
  right_flow = RightBoundFlow();
  for (int jVar = 0; jVar < _nbr_variables; jVar++)
  {
    flow_contribution[0][jVar] += left_flow[jVar];
    flow_contribution[_nbr_elements-1][jVar] -= right_flow[jVar];
  }

  //Calcul de U^(n+1) à partir de U^n
  double coeff = _delta_t/_delta_x;
  for (int elem_i = 0; elem_i < _nbr_elements; elem_i++)
  {
    for (int jVar = 0; jVar < _nbr_variables; jVar++)
    {
      _U[elem_i][jVar] -= coeff*flow_contribution[elem_i][jVar];
    }
  }
}

std::vector<double> Probleme1D::RightInterfaceFlow(int elem_i)
{
  std::vector<double> flow;

  switch (_choix_flux)
  {
    case left_elem_flow:
    return _flows[elem_i];
    break;

    case arithm_avg:
    flow.resize(_nbr_variables);
    for (int jVar = 0; jVar < _nbr_variables; jVar++)
    {
      flow[jVar] = (_flows[elem_i][jVar] + _flows[elem_i+1][jVar])/2.;
    }
    return flow;
    break;

    case harm_avg:
    flow.resize(_nbr_variables);
    for (int jVar = 0; jVar < _nbr_variables; jVar++)
    {
      flow[jVar] = 2./(1./_flows[elem_i][jVar] + 1./_flows[elem_i+1][jVar]);
    }
    return flow;
    break;

    default:
    std::cout << "Problème avec le choix de l'approximation des flux" << std::endl;
    exit(1);
  }
}

//Méthodes particulières
void TroisVariables::UpdateFlows()
{
  //retourne F_i+1/2

}

std::vector<double> TroisVariables::LeftBoundFlow()
{
  return {0,0,0};
}

std::vector<double> TroisVariables::RightBoundFlow()
{
  return {0,0,0};
}
