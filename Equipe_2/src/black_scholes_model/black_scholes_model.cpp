#include "black_scholes_model.hpp"
class black_scholes_model
{
private:
    /* data */
    double interest_rate;
    PnlVect* volatility;
    PnlVect* spots;
    double correlation;
public:
    black_scholes_model(/* args */);
    ~black_scholes_model();
};

black_scholes_model::black_scholes_model(double rate, PnlVect* vol, PnlVect* spots, double corr)
{
    interest_rate = rate;
    volatility = vol;
    spots = spots;
    correlation = corr;
}

black_scholes_model::~black_scholes_model()
{
    pnl_vect_free(&volatility);
    pnl_vect_free(&spots);
}
