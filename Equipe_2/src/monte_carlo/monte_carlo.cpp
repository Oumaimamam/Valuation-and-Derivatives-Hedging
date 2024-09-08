#include "monte_carlo.hpp"

PnlVect *MonteCarlo::getDates() const
{
  int size = this->fixing_dates_number + 1;
  PnlVect *vect_ti = pnl_vect_create(size);
  // Maturity
  double T = option->maturity;
  for (int k = 0; k < size; k++)
  {
    pnl_vect_set(vect_ti, k, k*T/this->fixing_dates_number);
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

MonteCarlo::~MonteCarlo() {
    // Libérer la mémoire si nécessaire
    if (option != nullptr) {
        delete option;
    }
    if (model != nullptr) {
        delete model;
    }
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


// Destructeur


// Méthode privée pour obtenir les dates
PnlVect* MonteCarlo::getDates() const {
    PnlVect* Dates = pnl_vect_create(fixing_dates_number+1);
    //pnl_vect_set(PnlVect ∗v, int i, double x) to set x in place i
    for (int k=0, k<fixing_dates_number+1, k++)
    {
        pnl_vect_set(Dates, k, k*option.getMaturity()/fixing_dates_number);
    }


}

// Méthode pour calculer le prix
double MonteCarlo::price(double t) {
    double r = model.interesetRate;
    int T = option.maturity;
    return exp(-r*T)*(1/sample_number)*sum();
}

double MonteCarlo::sum(){
    double res = 0.0;
    for(int j = 1, j<fixing_dates_number+1, j++)
    {
        pnlMat* vectSim = pnl_mat_create(1, fixing_dates_number+1);
        res += option.payoff(asset(getDates()),vectSim);
    }
    return res/sample_number ;
}
