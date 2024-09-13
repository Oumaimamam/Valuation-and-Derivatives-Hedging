#include <gtest/gtest.h>
#include "../Utils/convert.hpp"
#include "../monte_carlo/monte_carlo.hpp"
#include <iostream>

// test deltas option call t=0

TEST(MonteCarloTest, TestingDeltasCall)
{
    MonteCarlo *monte_carlo = convert_json_to_monte_carlo("../../data/call/call.json", "../../data/call/call_market.txt");

    int option_size = monte_carlo->option->option_size;
    PnlVect* deltas_vect = pnl_vect_create(option_size); // For storing deltas
    PnlVect* stddev_deltas_vect = pnl_vect_create(option_size);  // For storing standard deviations of deltas
    PnlVect* s_t = pnl_vect_new();  // Store the spot prices at t=0
    PnlMat* cots = pnl_mat_new();   // Store the market data for asset prices at t=0

    monte_carlo->get_cotations(0, cots, s_t);

    monte_carlo->delta(cots, deltas_vect, stddev_deltas_vect, 0.0);

    std::cout << "Deltas of the call option: ";
    for (int i = 0; i < deltas_vect->size; i++) {
        std::cout << pnl_vect_get(deltas_vect, i) << " ";
    }
    std::cout << std::endl;

    std::cout << "Expected deltas of the call option: 0.6314834286591736" << std::endl;

    std::cout << "Standard deviations of deltas: ";
    for (int i = 0; i < stddev_deltas_vect->size; i++) {
        std::cout << pnl_vect_get(stddev_deltas_vect, i) << " ";
    }
    std::cout << "Expected standard deviations of deltas: 0.00234292944876985" << std::endl; 
    std::cout << std::endl;

    delete monte_carlo;
    pnl_vect_free(&deltas_vect);
    pnl_vect_free(&stddev_deltas_vect);
    pnl_vect_free(&s_t);
    pnl_mat_free(&cots);
}

// #include <gtest/gtest.h>

// TEST(ExampleTest, HandlesTrueAssertions)
// {
//     EXPECT_TRUE(true);
// }

// TEST(ExampleTest, HandlesFalseAssertions)
// {
//     EXPECT_FALSE(false);
// }
