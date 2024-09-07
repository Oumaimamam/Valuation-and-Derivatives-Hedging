#include "option.hpp"

Option::Option(double T, double K, OptionType type, double D, PnlVect *coeff)
    : maturity(T),
      strike(K),
      type(type),
      option_size(D),
      payoff_coeffcients(coeff)
{
}

Option::Option(const Option &autre)
{
    this->maturity = autre.maturity;
    this->strike = autre.strike;
    this->type = autre.type;
    this->option_size = autre.option_size;
    this->payoff_coeffcients = pnl_vect_copy(autre.payoff_coeffcients);
}

Option::~Option()
{
    pnl_vect_free(&this->payoff_coeffcients);
}

Option *Option::GetOption(double T, double K, OptionType type, double D, PnlVect *coeff)
{
    switch (type)
    {
    case OptionType::Basket:
        Option *option = new OptionBasket(T, K, type, D, coeff);
        return option;

    case OptionType::Asian:
        Option *option = new OptionAsian(T, K, type, D, coeff);
        return option;

    case OptionType::Performance:
        Option *option = new OptionPerformance(T, K, type, D, coeff);
        return option;

    default:
        throw std::invalid_argument("Type d'option non valide : " + type);
    }
}
