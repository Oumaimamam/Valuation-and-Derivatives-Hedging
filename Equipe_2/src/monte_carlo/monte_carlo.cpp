#include "monte_carlo.hpp"
#include <cmath>

void MonteCarlo::get_all_dates(PnlVect *vect) const
{
    /*
  To use this function there are some steps to do :
  PnlVect* vect = pnl_vect_new();
  get_all_dates(vect);
  pnl_vect_free(&vect);

  */
    int size = this->fixing_dates_number + 1;
    pnl_vect_resize(vect, size);
    // Maturity
    double T = option->maturity;
    for (int k = 0; k < size; k++)
    {
        pnl_vect_set(vect, k, k * T / this->fixing_dates_number);
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
    // Libérer la mémoire si nécessaire
    if (option != nullptr)
    {
        delete option;
    }
    if (model != nullptr)
    {
        delete model;
    }
}

double MonteCarlo::price(double t)
{
    // calcul du prix à l'instat t = 0 et pour d = 1
    double v_0 = 0.0;

    PnlVect *dates = pnl_vect_new();
    get_all_dates(dates);


    PnlMat *matrix = pnl_mat_create(this->option->option_size, dates->size);

    for (int i = 1; i < this->sample_number + 1; i++)
    {
        this->model->asset(dates, matrix);
        v_0 += this->option->payOff(matrix);

    }

    double r = this->model->interest_rate;
    int T = this->option->maturity;

    // free :
    pnl_vect_free(&dates);
    pnl_mat_free(&matrix);

    return std::exp(-r * T) * (1 / this->sample_number) * v_0;
}

// Destructeur

// Méthode privée pour obtenir les dates
// PnlVect *MonteCarlo::getDates() const
// {
//     PnlVect *Dates = pnl_vect_create(fixing_dates_number + 1);
//     // pnl_vect_set(PnlVect ∗v, int i, double x) to set x in place i
//     for (int k = 0, k < fixing_dates_number + 1, k++)
//     {
//         pnl_vect_set(Dates, k, k * option.getMaturity() / fixing_dates_number);
//     }
// }

// Méthode pour calculer le prix
// double MonteCarlo::price(double t)
// {
//     double r = model.interesetRate;
//     int T = option.maturity;
//     return exp(-r * T) * (1 / sample_number) * sum();
// }

// double MonteCarlo::sum()
// {
//     double res = 0.0;
//     for (int j = 1, j < fixing_dates_number + 1, j++)
//     {
//         pnlMat *vectSim = pnl_mat_create(1, fixing_dates_number + 1);
//         res += option.payoff(asset(getDates()), vectSim);
//     }
//     return res / sample_number;
// }
