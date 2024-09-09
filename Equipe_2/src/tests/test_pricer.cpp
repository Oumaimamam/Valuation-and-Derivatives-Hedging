#include <gtest/gtest.h>
#include "../Utils/convert.hpp"
#include "../monte_carlo/monte_carlo.hpp"
#include <iostream>

TEST(MonteCarloTest, TestingPriceCall)
{
    MonteCarlo *monte_carlo = convert_json_to_monte_carlo("../../data/call/call.json");
    double price0 = monte_carlo->price(0);
    std::cout << "price0 = " << price0;
}

TEST(ExampleTest, HandlesFalseAssertions)
{
    EXPECT_FALSE(false);
}
