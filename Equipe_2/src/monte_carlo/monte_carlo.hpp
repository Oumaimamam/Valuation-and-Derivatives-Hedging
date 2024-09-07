#ifndef MONTE_CARLO_HPP
#define MONTE_CARLO_HPP

#include "../Option/option.hpp"
#include "../black_scholes_model/black_scholes_model.hpp"
#include "pnl/pnl_vector.h"
class MonteCarlo
{
public:
    Option *option;
    BlackScholesModel *model;
    int fixing_dates_number;
    int sample_number;
    void getDates(PnlVect *vect) const;

public:
    MonteCarlo(Option *option, BlackScholesModel *model, int N, int M);
    ~MonteCarlo();
    void price(double t);
};

#endif