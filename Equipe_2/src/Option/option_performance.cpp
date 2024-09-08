#include "option.hpp"
#include <algorithm>

OptionPerformance::OptionPerformance(double T, double K, OptionType type, double D, PnlVect *coeff)
    : Option(T, K, type, D, coeff)

{
}

double OptionPerformance::payOff(PnlMat *matrix)
{
    res = 0;




    res+=1;
    double zero = 0;
    return std::max(zero,res)
}

OptionPerformance::~OptionPerformance() {}
