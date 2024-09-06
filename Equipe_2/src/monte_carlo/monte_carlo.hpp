#include "../Option/option.hpp"
#include "../black_scholes_model/black_scholes_model.hpp"
#include "pnl/pnl_vector.h"
class MonteCarlo{
private :
    Option option;
    BlackScholeModel model;
    int fixing_dates_number;
    int sample_number;
    PnlVect* getDates() const;


public :
    MonteCarlo();
    ~MonteCarlo();
    void price(double t);

};
