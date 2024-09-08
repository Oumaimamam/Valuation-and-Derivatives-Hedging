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
    PnlVect *volatility; //pnl_vect_free
    PnlVect *spots; //pnl_vect_free
    double correlation;
    PnlMat *mat_asset ;  //pnl_mat_free
    int model_size;
    PnlRng *rng;  //pnl_rng_free(&rng)

public:
    BlackScholesModel();
    BlackScholesModel(double rate, PnlVect *vol, PnlVect *spots, double corr);
    ~BlackScholesModel();
    void asset(PnlVect *Dates);
    void get_matrix_Cholesky_corr(PnlMat *matrix);
};
#endif