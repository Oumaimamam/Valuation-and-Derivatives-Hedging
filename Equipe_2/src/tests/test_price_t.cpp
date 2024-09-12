#include <iostream>
#include "pnl/pnl_cdf.h"
#include <cmath>
#include <fstream>
#include "pnl/pnl_vector.h"
#include "../monte_carlo/monte_carlo.hpp"
#include "../Utils/convert.hpp"
#include "../Utils/utils.hpp"

double calcul_vt_expected(double t, double s_t, double k, double r, double T, double sigma)
{
    double d_1 = (std::log(s_t / k) + (r + pow(sigma, 2)) * (T - t)) / (sigma * sqrt(T - t));
    double d_2 = d_1 - sigma * sqrt(T - t);

    double v_t = s_t * pnl_cdfnor(d_1) - k * exp(-r * (T - t)) * pnl_cdfnor(d_2);

    return v_t;
}

void get_cotations(double t, PnlMat *past, PnlMat *market_data, MonteCarlo *monte_carlo)
{
    /*
    To use this function there are some steps to do :

    Step 1 :
    PnlVect* s_t = pnl_vect_new();
    PnlMat* cots = pnl_mat_new();

    get_cotations(t , cots , s_t);

    pnl_vect_free(&s_t);
    pnl_mat_free(&cots);
  */

    if (market_data == NULL)
    {
        throw std::invalid_argument("argument data.txt non fourni");
        exit(1);
    }

    int H = monte_carlo->model->hedging_dates_number;
    int N = monte_carlo->fixing_dates_number;
    double T = monte_carlo->option->maturity;
    int D = monte_carlo->option->option_size;

    int i = compute_last_index(t, T, N);

    pnl_mat_resize(past, i + 2, D);

    PnlVect *col = pnl_vect_create(D);

    for (int j = 0; j < i + 1; j++)
    {
        pnl_mat_get_row(col, market_data, j * H / N);
        pnl_mat_set_row(past, col, j);
    }

    // s_t
    pnl_mat_get_row(col, market_data, t * H / T);
    pnl_mat_set_row(past, col, i + 1);

    // free
    pnl_vect_free(&col);
}

int main()
{
    MonteCarlo *monte_carlo = convert_json_to_monte_carlo("../../data/call/call.json");

    double K = monte_carlo->option->strike;
    double r = monte_carlo->model->interest_rate;
    double T = monte_carlo->option->maturity;
    double sigma = monte_carlo->model->volatility->array[0];

    double H = monte_carlo->model->hedging_dates_number;
    int i = 10;
    double t = i * T / H;

    PnlMat *data = pnl_mat_create_from_file("../../data/call/call_market.txt");

    PnlMat *past = pnl_mat_new();

    get_cotations(t, past, data, monte_carlo);

    int D = monte_carlo->option->option_size;
    PnlVect *vect_st = pnl_vect_create(D);
    pnl_mat_get_row(vect_st, past, past->m - 1);

    double s_t = vect_st->array[0]; // s_t , D = 1

    // v_t expected :

    double v_t_expected = calcul_vt_expected(t, s_t, K, r, T, sigma);

    // calcul du v_t (pricer) :
    PnlMat *matrix = pnl_mat_new();
    double v_t;
    double v_t_std;
    monte_carlo->price(t, v_t, v_t_std, past, matrix);

    std::cout << "v_t  =" << v_t << std::endl;
    std::cout << "v_t_expected  =" << v_t_expected << std::endl;

    pnl_mat_free(&data);
    pnl_mat_free(&past);
    pnl_vect_free(&vect_st);
    pnl_mat_free(&matrix);
    delete monte_carlo;

    return 0;
}