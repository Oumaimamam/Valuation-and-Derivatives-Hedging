#include "monte_carlo.hpp"

PnlVect *MonteCarlo::getDates() const
{
  int size = this->fixing_dates_number + 1;
  PnlVect *vect_ti = pnl_vect_create(size);
  // Maturity
  double T;
  for (int k = 0; k < size; k++)
  {
    vect_ti->array[k] = (k * T) / (this->fixing_dates_number);
  }

  return vect_ti;
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
