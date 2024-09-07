#include <gtest/gtest.h>
#include "../monte_carlo/monte_carlo.hpp"
#include <iostream>
#include <fstream>
#include <nlohmann/json.hpp>
#include "../Utils/convert.hpp"

TEST(MonteCarloTest, TestingGetDatesDataAsianJSON)
{
    std::ifstream file("../../data/asian/asian.json");
    if (!file.is_open())
    {
        std::cerr << "Error opening file" << std::endl;
        exit(1);
    }
    nlohmann::json json = nlohmann::json::parse(file);

    MonteCarlo *monte_carlo = convert_json_to_monte_carlo(json);
    PnlVect *list_ti = pnl_vect_new();
    monte_carlo->getDates(list_ti);

    EXPECT_TRUE(list_ti->size == 25);
    EXPECT_TRUE(list_ti->array[0] == 0.0);
    EXPECT_TRUE(list_ti->array[24] == 1.5);

    file.close();
    pnl_vect_free(&list_ti);
}