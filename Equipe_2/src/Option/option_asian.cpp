#include "option.hpp"
#include <algorithm>

OptionAsian::OptionAsian(double T, double K, OptionType type, double D, PnlVect *coeff)
    : Option(T, K, type, D, coeff)

{
}

double OptionAsian::payOff(PnlMat *matrix)
{
    double res = 0;
    int rows = matrix->m; // D
    int cols = matrix->n; // N+1

    for (int d = 0; d < rows; d++)
    {
        double lamda_d = GET(payoff_coeffcients, d);
        double res2 = 0;
        for (int j = 0; j < cols; j++)
        {
            res2 += pnl_mat_get(matrix, d, j);
        }

        res += lamda_d * res2;
    }

    res = res / (cols) - this->strike;
    double zero = 0.0;
    return std::max(zero, res);
}
OptionAsian::~OptionAsian() {}
