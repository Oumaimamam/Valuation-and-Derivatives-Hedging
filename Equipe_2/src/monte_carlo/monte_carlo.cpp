#include "monte_carlo.hpp"
#include <cmath>
#include "../Utils/convert.hpp"
#include "../Utils/utils.hpp"

void MonteCarlo::get_all_dates(PnlVect *vect, double t, int i) const
{
    /*
  To use this function there are some steps to do :
  PnlVect* vect = pnl_vect_new();
  get_all_dates(vect);
  pnl_vect_free(&vect);

  */
    int size = this->fixing_dates_number + 1 - i;
    pnl_vect_resize(vect, size);
    double T = option->maturity;

    pnl_vect_set(vect, 0, 0);
    for (int k = 1; k < size; k++)
    {
        double u_k = ((k + i) * T / this->fixing_dates_number) - t;
        pnl_vect_set(vect, k, u_k);
    }
}

MonteCarlo::MonteCarlo(Option *option, BlackScholesModel *model, int N, int M, PnlMat *data, double h)
    : option(option),
      model(model),
      fixing_dates_number(N),
      sample_number(M),
      market_data(data),
      fd_step(h)
{
}

MonteCarlo::~MonteCarlo()
{
    // Libérer la mémoire si nécessaire
    if (option != nullptr)
    {
        delete option;
    }
    if (model != nullptr)
    {
        delete model;
    }

    if (market_data != nullptr)
    {
        pnl_mat_free(&market_data);
    }
}

void MonteCarlo::price(double t, double &price, double &price_std)
{

    int D = this->option->option_size;
    double r = this->model->interest_rate;
    int M = this->sample_number;
    int N = this->fixing_dates_number;
    double T = this->option->maturity;


    double v_t = 0.0;
    double price_std_dev = 0.0;


    PnlMat *matrix = pnl_mat_create(D, N + 1);



    for (int i = 0; i < M; i++)
    {
        get_matrix_of_sim(t , matrix);
        double phi_j = this->option->payOff(matrix);
        v_t += phi_j;
        price_std_dev += pow(phi_j, 2);
    }


    double inv_M = 1.0 / (double)M;

    price = std::exp(-r * (T - t)) * inv_M * v_t;

    price_std = sqrt(exp(-2 * r * T) * (inv_M * price_std_dev - pow(inv_M * v_t, 2))) / sqrt(M);

    pnl_mat_free(&matrix);
}

void MonteCarlo::get_cotations(double t, PnlMat *cots, PnlVect *s_t)
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

    int H = this->model->hedging_dates_number;
    int N = this->fixing_dates_number;
    double T = this->option->maturity;
    int D = this->option->option_size;

    int i = compute_last_index(t, T, N);

    pnl_mat_resize(cots, D, i + 1);

    PnlVect *col = pnl_vect_create(D);

    for (int j = 0; j < i + 1; j++)
    {
        pnl_mat_get_row(col, this->market_data, j * H / N);
        pnl_mat_set_col(cots, col, j);
    }

    int index_t = t * H / T;
    pnl_mat_get_row(s_t, this->market_data, index_t);

    pnl_vect_free(&col);
}



void MonteCarlo::get_matrix_of_sim(double t , PnlMat *matrix)
{

    double T = this->option->maturity;
    int N = this->fixing_dates_number;
    int index = compute_last_index(t, T, N);
    int D = this->option->option_size;

    PnlVect *dates = pnl_vect_new();
    get_all_dates(dates, t, index);

    PnlMat *matrix_sim = pnl_mat_create(D, N - index);


    if (t == 0.0)
    {
        this->model->asset(t, dates, matrix_sim);
        pnl_mat_set_col(matrix , this->model->spots , 0);
        pnl_mat_set_subblock(matrix, matrix_sim, 0, index + 1);
        
        
    }
    else
    {
        PnlVect *s_t = pnl_vect_new();
        PnlMat *cots = pnl_mat_new();
        get_cotations(t, cots, s_t);

        PnlVect *col = pnl_vect_create(D);

        this->model->asset(t, dates, matrix_sim);

        for (int j = 0; j < N - index; j++)
        {
            pnl_mat_get_col(col, matrix_sim, j);
            pnl_vect_mult_vect_term(col, s_t);
            pnl_mat_set_col(matrix_sim, col, j);
        }

        pnl_mat_set_subblock(matrix, cots, 0, 0);
        pnl_mat_set_subblock(matrix, matrix_sim, 0, index + 1);

        pnl_vect_free(&col);
        pnl_vect_free(&s_t);
        pnl_mat_free(&cots);
    }

    pnl_mat_free(&matrix_sim);
    pnl_vect_free(&dates);

}




// PnlVect* MonteCarlo::delta(PnlVect* delta_vect)
// {
//     int D = this->option->option_size; 
//     int N = this->fixing_dates_number; 
//     double T = this->option->maturity;
//     int M = this->sample_number; 
//     double h = this->fd_step;
//     int i = compute_last_index(t, T, N);
//     pnl_vect_resize(deltas, D);

//     PnlMat *cots = pnl_mat_new(); 
//     PnlVect *s_t = pnl_vect_new(); 
//     PnlMat *mat_plus_h = pnl_mat_new(); 
//     PnlMat *mat_minus_h = pnl_mat_new(); 
//     PnlMat *matrix_sim = pnl_mat_create(D, this->fixing_dates_number - index);
//     get_cotations(t, cots, s_t);

//     for (int d = 0; d < D; d++)
//     {
//         double delta_sum = 0.0;
//         for (int j = 0; j < M; j++)
//         {
//             pnl_vect_clone(s_t, this->model->spots);
//             pnl_vect_set(s_t, d, pnl_vect_get(s_t, d) * (1 + h));
//             this->model->asset(s_t, cots, mat_plus_h);

//             pnl_vect_set(s_t, d, pnl_vect_get(s_t, d) * (1 - 2 * h));
//             this->model->asset(s_t, cots, mat_minus_h);

//             double payoff_plus = this->option->payOff(mat_plus_h);
//             double payoff_minus = this->option->payOff(mat_minus_h);

            
//             delta_sum += (payoff_plus - payoff_minus) / (2 * h);
//         }

        
//         double delta_d = delta_sum / M;
//         pnl_vect_set(deltas, d, delta_d);
//     }
//     pnl_mat_free(&cots);
//     pnl_vect_free(&s_t);
//     pnl_mat_free(&mat_plus_h);
//     pnl_mat_free(&mat_minus_h);
// }
