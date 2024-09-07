#include "black_scholes_model.hpp"

BlackScholesModel::BlackScholesModel(double rate, PnlVect *vol, PnlVect *spots, double corr)
    : interest_rate(rate),
      volatility(vol),
      spots(spots),
      correlation(corr)
{
}

BlackScholesModel::~BlackScholesModel()
{
    pnl_vect_free(&this->volatility);
    pnl_vect_free(&this->spots);
}
