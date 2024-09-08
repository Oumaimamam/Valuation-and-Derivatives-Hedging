#include "option.hpp"

OptionBasket::OptionBasket(double T, double K, OptionType type, double D, PnlVect *coeff)
    : Option(T, K, type, D, coeff)

{
}

double OptionBasket::payOff(PnlMat *matrix)
{
    double res = 0;
    int rows,col = matrix->m;
    int cols = matrix ->n;

    for(int i =0; i<rows;i++)
    {
        double lamda_d = GET(payoff_coeffcients,i);
        double S_t_d = pnl_mat_get(matrix,i,cols);
        
        res += lamda_d*S_t_d;
    }
    res -= this->strike;
    
    double zero = 0.0;
    return std::max(zero,res)
}

OptionBasket::~OptionBasket() {}

