#ifndef MONTE_CARLO_HPP
#define MONTE_CARLO_HPP

#include "../Option/option.hpp"
#include "../black_scholes_model/black_scholes_model.hpp"
#include "pnl/pnl_vector.h"
#include "../pcpd_helper.hpp"

class MonteCarlo
{
public:
    Option *option;
    BlackScholesModel *model;
    int fixing_dates_number;
    int sample_number;
    PnlMat* market_data ; 
    // PricingResults res_price ;
    // HedgingResults  res_hedge ; 

    void get_all_dates(PnlVect *vect , double t , int i ) const;

public:
    MonteCarlo(Option *option, BlackScholesModel *model, int N, int M , PnlMat* data );
    ~MonteCarlo();
    // calculer le price
    double price(double t);
    void get_cotations(double t, PnlMat *cots ,  PnlVect* s_t);
};

#endif