#include <gtest/gtest.h>
#include "../monte_carlo/monte_carlo.hpp"
#include <iostream>
#include <fstream>
#include <nlohmann/json.hpp>
#include "../Utils/convert.hpp"

TEST(MonteCarloTest, TestingGetDatesDataAsianJSON)
{

    MonteCarlo *monte_carlo = convert_json_to_monte_carlo("../../data/asian/asian.json");
    PnlVect *list_ti = pnl_vect_new();
    monte_carlo->get_all_dates(list_ti, 0, 0);

    EXPECT_TRUE(list_ti->size == 25);
    EXPECT_TRUE(list_ti->array[0] == 0.0);
    EXPECT_TRUE(list_ti->array[2] == 0.125);
    EXPECT_TRUE(list_ti->array[24] == 1.5);

    pnl_vect_free(&list_ti);
    delete monte_carlo;
}