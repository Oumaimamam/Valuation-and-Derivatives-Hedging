#include "option.hpp"
#include <algorithm>
OptionBasket::OptionBasket(double T, double K, OptionType type, double D, PnlVect *coeff)
    : Option(T, K, type, D, coeff)

{
}

double OptionBasket::payOff(PnlMat *matrix)
{
    double res = 0;
    int rows = matrix->m;
    int cols = matrix ->n;

    for(int d =0; d<rows;d++)
    {
        double lamda_d = GET(payoff_coeffcients,d);
        double S_t_d = pnl_mat_get(matrix,d,cols-1);
        
        res += lamda_d*S_t_d;
    }
    res -= this->strike;
    
    double zero = 0.0;
    return std::max(zero,res);
}

OptionBasket::~OptionBasket() 
{

}

