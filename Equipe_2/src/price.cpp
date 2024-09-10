#include <iostream>
#include "monte_carlo/monte_carlo.hpp"
#include "Utils/convert.hpp"
#include "monte_carlo/monte_carlo.hpp"
#include "pnl/pnl_vector.h"

int main(int argc, char *argv[])
{
    if(argc != 2) {
        std::cout << "" << std::endl ; 
    }

    MonteCarlo* monte_carlo = convert_json_to_monte_carlo(argv[1]);

    double prix = monte_carlo->price(0);
    double prix_std_dev ; 
    PnlVect* delta ;
    PnlVect* delta_std_dev ;



    PricingResults res(prix, prix_std_dev, delta, delta_std_dev);
    std::cout << res << std::endl;

    delete monte_carlo;

    return 0;
}