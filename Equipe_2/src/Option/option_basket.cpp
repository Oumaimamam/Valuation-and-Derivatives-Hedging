#include "option.hpp"

OptionBasket::OptionBasket(double T, double K, OptionType type, double D, PnlVect *coeff)
    : Option(T, K, type, D, coeff)

{
}

void OptionBasket::payOff(PnlMat *matrix)
{
}

OptionBasket::~OptionBasket() {}
