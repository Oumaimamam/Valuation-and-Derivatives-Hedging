#ifndef BLACK_SCHOLES_MODEL_HPP
#define BLACK_SCHOLES_MODEL_HPP
#include "pnl/pnl_matvect.h"
#include "pnl/pnl_vector.h"
#include "pnl/pnl_matrix.h"
#include "pnl/pnl_random.h"
class BlackScholesModel
{
public:
    double interest_rate;
    PnlVect *volatility; // pnl_vect_free
    PnlVect *spots;      // pnl_vect_free
    double correlation;
    int model_size;
    int hedging_dates_number;
    void get_all_dates(PnlVect *vect, double t, int i) const;

public:
    BlackScholesModel();
    BlackScholesModel(double rate, PnlVect *vol, PnlVect *spots, double corr, double H);
    ~BlackScholesModel();
    // void asset(double t, PnlVect *Dates, PnlMat *mat_asset);
    void asset(const PnlMat *past, double t, PnlMat *path, PnlRng *rng, double maturity);
    void get_matrix_Cholesky_corr(PnlMat *matrix);
};
#endif