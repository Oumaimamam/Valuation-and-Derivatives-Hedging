#include <gtest/gtest.h>
#include "../Utils/convert.hpp"
#include "../monte_carlo/monte_carlo.hpp"
#include <iostream>

// test deltas option call

TEST(MonteCarloTest, TestingDeltasCall)
{
    MonteCarlo *monte_carlo = convert_json_to_monte_carlo("../../data/call/call.json", "../../data/call/call_market.txt");
    PnlVect* deltas_vect = pnl_vect_new();
    PnlVect* s_t = pnl_vect_new();
    PnlMat* cots = pnl_mat_new();
    monte_carlo->get_cotations(0, cots, s_t);
    monte_carlo->delta(deltas_vect, 0, cots, s_t);
    std::cout << "deltas of call option : "; 
    for (int i = 0; i < deltas_vect -> size; i++){
        std::cout << pnl_vect_get(deltas_vect, i) << " ";
    }
    std::cout << std::endl;
    
    delete monte_carlo;
    pnl_vect_free(&deltas_vect);
    pnl_vect_free(&s_t);
    pnl_mat_free(&cots);

}