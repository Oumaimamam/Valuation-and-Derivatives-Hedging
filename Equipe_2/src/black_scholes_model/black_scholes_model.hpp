#ifndef BLACK_SCHOLES_MODEL_HPP
#define BLACK_SCHOLES_MODEL_HPP
#include "pnl/pnl_matvect.h"
#include "pnl/pnl_vector.h"
class BlackScholesModel
{
public:
    double interest_rate;
    PnlVect *volatility;
    PnlVect *spots;
    double correlation;
    int model_size;

public:
    BlackScholesModel();
    BlackScholesModel(double rate, PnlVect *vol, PnlVect *spots, double corr);
    ~BlackScholesModel();
    void asset();
    void get_matrix_Cholesky_corr(PnlMat *matrix);
};
#endif