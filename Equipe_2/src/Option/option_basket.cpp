#include "option.hpp"
#include <algorithm>
OptionBasket::OptionBasket(double T, double K, OptionType type, double D, PnlVect *coeff)
    : Option(T, K, type, D, coeff)

{
}

double OptionBasket::payOff(PnlMat *matrix)
{
    double res = 0.0;
    int D = this->option_size;
    
    PnlVect* last_col = pnl_vect_create(D);
    pnl_mat_get_col(last_col , matrix , matrix->n -1);
    res = pnl_vect_scalar_prod(last_col , this->payoff_coeffcients) - this->strike;
    pnl_vect_free(&last_col);
    return std::max(0.0, res);

}

OptionBasket::~OptionBasket()
{
}
