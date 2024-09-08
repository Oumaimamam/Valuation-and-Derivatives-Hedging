#include "../Option/option.hpp"
#include "../black_scholes_model/black_scholes_model.hpp"
#include "pnl/pnl_vector.h"
class MonteCarlo{
private :
    Option *option;
    BlackScholesModel *model;
    int fixing_dates_number;
    int sample_number;
    void get_all_dates(PnlVect *vect) const;


public :
    MonteCarlo();
    ~MonteCarlo();
    // calculer le price
    double price(double t);
};
