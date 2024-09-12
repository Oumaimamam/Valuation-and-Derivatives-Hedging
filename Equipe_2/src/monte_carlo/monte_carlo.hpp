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
    double fd_step;
    PnlRng *rng; 
    int hedging_date_number;
    // PricingResults res_price ;
    // HedgingResults  res_hedge ; 

    // void get_all_dates(PnlVect *vect, double t, int i) const;

public:
    MonteCarlo(Option *option, BlackScholesModel *model, int N, int M, double h);
    ~MonteCarlo();
    // calculer le price
    void price(double &price, double &price_std);
    void price(double t , double &price, double &price_std , const PnlMat* Past, PnlMat *matrix);
    void calculPAndL(const PnlMat *past);

    void get_cotations(double t, PnlMat *cots, PnlMat *market_data);
    // void get_matrix_of_sim(double t , PnlMat *matrix);
    // void couverture

    // calculer les deltas
    void MonteCarlo::get_delta(PnlVect* Vect);

    void getCouvPrtf(PnlMat* Past, PnlMat* market_data, double &p_and_l, double t);
};

#endif