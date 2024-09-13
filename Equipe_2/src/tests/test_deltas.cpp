#include <gtest/gtest.h>
#include "../Utils/convert.hpp"
#include "../monte_carlo/monte_carlo.hpp"
#include "../Utils/compute_last_index.hpp"
#include <iostream>

// test deltas option call t=0
// void get_cotations(double t, PnlMat *past, PnlMat *market_data, MonteCarlo *monte_carlo)
// {
//     /*
//     To use this function there are some steps to do :

//     Step 1 :
//     PnlVect* s_t = pnl_vect_new();
//     PnlMat* cots = pnl_mat_new();

//     get_cotations(t , cots , s_t);

//     pnl_vect_free(&s_t);
//     pnl_mat_free(&cots);
//   */

//     int H = monte_carlo->model->hedging_dates_number;
//     int N = monte_carlo->fixing_dates_number;
//     double T = monte_carlo->option->maturity;
//     int D = monte_carlo->option->option_size;

//     int i = compute_last_index(t, T, N);

//     pnl_mat_resize(past, i + 2, D);

//     PnlVect *col = pnl_vect_create(D);

//     for (int j = 0; j < i + 1; j++)
//     {
//         pnl_mat_get_row(col, market_data, j * H / N);
//         pnl_mat_set_row(past, col, j);
//     }

//     // s_t
//     pnl_mat_get_row(col, market_data, t * H / T);
//     pnl_mat_set_row(past, col, i + 1);

//     // free
//     pnl_vect_free(&col);
// }

int main()
{
    MonteCarlo *monte_carlo = convert_json_to_monte_carlo("../../data/call/call.json");

    int option_size = monte_carlo->option->option_size;
    PnlVect* deltas_vect = pnl_vect_create(option_size);
    PnlVect* stddev_deltas_vect = pnl_vect_create(option_size);

    monte_carlo->delta(deltas_vect, stddev_deltas_vect);

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
    return 1;

}
