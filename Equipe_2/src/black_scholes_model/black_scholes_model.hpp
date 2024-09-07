#ifndef BLACK_SCHOLES_MODEL_HPP
#define BLACK_SCHOLES_MODEL_HPP
#include "pnl/pnl_matvect.h"
#include "pnl/pnl_vector.h"
class BlackScholesModel
{
private:
    double interest_rate;
    PnlVect *volatility;
    PnlVect *spots;
    double correlation;

public:
    BlackScholesModel();
    BlackScholesModel(double rate, PnlVect *vol, PnlVect *spots, double corr);
    ~BlackScholesModel();
    void asset();
};
#endif