#include <gtest/gtest.h>
#include "../monte_carlo/monte_carlo.hpp"
#include "../Utils/convert.hpp"

TEST(BlackScholesModel, TestAssetCall)
{
    MonteCarlo* monte_carlo = convert_json_to_monte_carlo("../../data/call/call.json");

    BlackScholesModel* model = monte_carlo->model;

}


