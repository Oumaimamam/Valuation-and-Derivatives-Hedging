#include "monte_carlo.hpp"

void MonteCarlo::getDates(PnlVect *vect) const
{
  int size = this->fixing_dates_number + 1;
  pnl_vect_resize(vect, size);
  // Maturity
  double T = option->maturity;
  for (int k = 0; k < size; k++)
  {
    vect->array[k] = (k * T) / (this->fixing_dates_number);
  }
}

MonteCarlo::MonteCarlo(Option *option, BlackScholesModel *model, int N, int M)
    : option(option),
      model(model),
      fixing_dates_number(N),
      sample_number(M)
{
}

MonteCarlo::~MonteCarlo()
{
  delete model;
  delete option;
}

void MonteCarlo::price(double t)
{
  // calcul du prix à l'instat t = 0 et pour d = 1
  double v_0 = 0.0;

  PnlVect *vect_phi_j;

  for (int i = 0; i < this->sample_number; i++)
  {
  }
}
