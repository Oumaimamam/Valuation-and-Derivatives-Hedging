#include "option.hpp"
#include <algorithm>

OptionPerformance::OptionPerformance(double T, double K, OptionType type, double D, PnlVect *coeff)
    : Option(T, K, type, D, coeff)

{
}

double OptionPerformance::payOff(PnlMat *matrix)
{
    int rows = matrix->m; // N+1
    int cols = matrix->n; // D
    PnlVect *S1_ti = pnl_vect_create(cols);
    double sum_prec = 0.0;
    double sum_curr = 0.0;
    double res = 0.0;
    pnl_mat_get_row(S1_ti, matrix, 0); //
    sum_prec = pnl_vect_scalar_prod(this->payoff_coeffcients, S1_ti);
    for (int i = 1; i < rows; i++)
    {
        pnl_mat_get_row(S1_ti, matrix, i); //
        sum_curr = pnl_vect_scalar_prod(this->payoff_coeffcients, S1_ti);
        res += std::max(sum_curr / sum_prec - 1.0, 0.0);
        sum_prec = sum_curr;
    }

    // free
    pnl_vect_free(&S1_ti);

    return 1 + res;
}

OptionPerformance::~OptionPerformance() {}

// double res = 0;
// int rows = matrix->m; // D
// int cols = matrix->n; // N+1

// for (int j = 1; j < cols; j++)
// {
//     double sum1 = 0;
//     double sum2 = 0;
//     double lamda_d;
//     for (int d = 0; d < rows; d++)
//     {
//         lamda_d = GET(payoff_coeffcients, d);
//         sum1 += lamda_d * pnl_mat_get(matrix, d, j);
//         sum2 += lamda_d * pnl_mat_get(matrix, d, j - 1);
//     }

//     res += std::max(sum1 / sum2 - 1, 0.0);
// }