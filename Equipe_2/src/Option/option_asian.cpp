#include "option.hpp"

OptionAsian::OptionAsian(double T, double K, OptionType type, double D, PnlVect *coeff)
    : Option(T, K, type, D, coeff)

{
}

void OptionAsian::payOff(PnlMat *matrix)
{
}

OptionAsian::~OptionAsian() {}
