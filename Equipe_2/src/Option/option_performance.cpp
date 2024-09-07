#include "option.hpp"

OptionPerformance::OptionPerformance(double T, double K, OptionType type, double D, PnlVect *coeff)
    : Option(T, K, type, D, coeff)

{
}

void OptionPerformance::payOff(PnlMat *matrix)
{
}

OptionPerformance::~OptionPerformance() {}
