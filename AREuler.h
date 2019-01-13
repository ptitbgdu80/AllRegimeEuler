#include <vector>
#include <math.h>
#include <iostream>

enum : int {left_elem_flow, arithm_avg, harm_avg};

class Probleme1D
{

protected:
  std::vector<std::vector<double>> _U, _flows;
  double _nbr_elements, _nbr_variables, _delta_x, _delta_t;
  int _choix_flux;

public:
  Probleme1D();
  ~Probleme1D();
  std::vector<double> RightInterfaceFlow(int elem_i);
  virtual std::vector<double> LeftBoundFlow() = 0;
  virtual std::vector<double> RightBoundFlow() = 0;
  virtual void UpdateFlows() = 0;
  void ClassicGodunovIteration();

};


class TroisVariables : public Probleme1D
{

public:
  TroisVariables(int nbr_elements, double delta_t, int choix_flux);
  std::vector<double> LeftBoundFlow();
  std::vector<double> RightBoundFlow();
  void UpdateFlows();
};
