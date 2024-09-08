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
    PnlMat *mat_asset ;

public:
    BlackScholesModel();
    BlackScholesModel(double rate, PnlVect *vol, PnlVect *spots, double corr);
    ~BlackScholesModel();
    void asset(PnlVect *Dates);
};
#endif