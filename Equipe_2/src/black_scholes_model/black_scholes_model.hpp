#ifndef BLACK_SCHOLES_MODEL_HPP
#define BLACK_SCHOLES_MODEL_HPP
#include "pnl/pnl_matvect.h"
#include "pnl/pnl_vector.h"
class black_scholes_model
{
    private:
        double interest_rate;
        PnlVect* volatility;
        PnlVect* spots;
        double correlation;

    public:
        black_scholes_model();
        ~black_scholes_model();
        void asset();

};
#endif