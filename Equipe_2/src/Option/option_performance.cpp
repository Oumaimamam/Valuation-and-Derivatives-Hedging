#include "option.hpp"
#include <algorithm>

OptionPerformance::OptionPerformance(double T, double K, OptionType type, double D, PnlVect *coeff)
    : Option(T, K, type, D, coeff)

{
}

double OptionPerformance::payOff(PnlMat *matrix)
{
    double res = 0;
    int rows = matrix->m; // D
    int cols = matrix->n; // N+1

    for (int j = 1; j < cols; j++)
    {
        double sum1 = 0;
        double sum2 = 0;
        double lamda_d;
        for (int d = 0; d < rows; d++)
        {
            lamda_d = GET(payoff_coeffcients, d);
            sum1 += lamda_d * pnl_mat_get(matrix, d, j);
            sum2 += lamda_d * pnl_mat_get(matrix, d, j - 1);
        }

        res += std::max(sum1 / sum2 - 1, 0.0);
    }

    return 1 + res;
}

OptionPerformance::~OptionPerformance() {}
