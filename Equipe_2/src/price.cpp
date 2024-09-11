#include <iostream>
#include "monte_carlo/monte_carlo.hpp"
#include "Utils/convert.hpp"
#include "monte_carlo/monte_carlo.hpp"
#include "pnl/pnl_vector.h"

int main(int argc, char *argv[])
{
    // if (argc != 2)
    // {
    //     std::cout << "" << std::endl;
    // }

    MonteCarlo *monte_carlo = convert_json_to_monte_carlo("../../data/basket/basket_40d/basket_40d.json", "../../data/basket/basket_40d/basket_40d_market.txt");

    double prix;
    double prix_std_dev;
    // PnlVect *delta;
    // PnlVect *delta_std_dev;

    monte_carlo->price(0.105233, prix, prix_std_dev);

    std::cout << "price =" << prix << std::endl;
    std::cout << "expected_t =" << 9.1280832137086 << std::endl;


    // PricingResults res(prix, prix_std_dev, delta, delta_std_dev);
    // std::cout << res << std::endl;

    delete monte_carlo;
    // pnl_vect_free(&delta);
    // pnl_vect_free(&delta_std_dev);

    return 0;
}