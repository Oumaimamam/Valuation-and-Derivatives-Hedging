#include <gtest/gtest.h>
#include "../Utils/convert.hpp"
#include "../monte_carlo/monte_carlo.hpp"
#include <iostream>

TEST(MonteCarloTest, TestingPriceCall)
{
    MonteCarlo *monte_carlo = convert_json_to_monte_carlo("../../data/call/call.json" , "../../data/call/call_market.txt");
    double price0 = monte_carlo->price(0.006703);
    std::cout << "price0 = " << price0 << std::endl;
}

TEST(MonteCarloTest, TestingPriceAsian)
{
    MonteCarlo *monte_carlo = convert_json_to_monte_carlo("../../data/asian/asian.json" , "../../data/asian/asian_market.txt");
    double price0 = monte_carlo->price(0.105233);
    std::cout << "price0 = " << price0 << std::endl;
}


// TEST(MonteCarloTest, TestingPriceAsian)
// {
//     MonteCarlo *monte_carlo = convert_json_to_monte_carlo("../../data/asian/asian.json" , "../../data/asian/asian_market.txt");
//     double price0 = monte_carlo->price(0.105233);
//     std::cout << "price0 = " << price0 << std::endl;
// }

