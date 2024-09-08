#include "option.hpp"
#include <algorithm>

OptionAsian::OptionAsian(double T, double K, OptionType type, double D, PnlVect *coeff)
    : Option(T, K, type, D, coeff)

{
}


double OptionAsian::payOff(PnlMat *matrix)
{
    double res = 0;
    int rows,col = matrix->m;
    int cols = matrix ->n;

    for(int i =0; i<rows;i++)
    {
        double lamda_d = GET(payoff_coeffcients,i);
        double res2=0;
        for (int j=0; j<cols+1; j++)
        {
            res2 += pnl_mat_get(matrix,i,j);
        }


        
        res += lamda_d*res2;
    }

    res = res/(rows) - this->strike;
    double zero = 0;
    return std::max(zero,res);
}
OptionAsian::~OptionAsian() {}
