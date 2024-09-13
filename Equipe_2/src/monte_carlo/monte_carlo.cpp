#include "monte_carlo.hpp"
#include <cmath>
#include "../Utils/convert.hpp"
// #include "../Utils/construct_append_mat.hpp"
#include "../Utils/utils.hpp"

MonteCarlo::MonteCarlo(Option *option, BlackScholesModel *model, int N, int M, double h)
    : option(option),
      model(model),
      fixing_dates_number(N),
      sample_number(M),
      fd_step(h)
{
    rng = pnl_rng_create(PNL_RNG_MERSENNE);
    pnl_rng_sseed(rng, time(NULL));
}

MonteCarlo::~MonteCarlo()
{
    delete option;
    delete model;
    pnl_rng_free(&rng);
}

void MonteCarlo::price(double &price, double &price_std)
{

    int D = this->option->option_size;
    double r = this->model->interest_rate;
    int M = this->sample_number;
    int N = this->fixing_dates_number;
    double T = this->option->maturity;

    double v_t = 0.0;
    double price_std_dev = 0.0;

    PnlMat *matrix = pnl_mat_create(N + 1, D);

    for (int i = 0; i < M; i++)
    {
        this->model->asset(matrix, this->rng);
        double phi_j = this->option->payOff(matrix);
        v_t += phi_j;
        price_std_dev += pow(phi_j, 2);
    }

    double inv_M = 1.0 / (double)M;

    price = std::exp(-r * (T)) * inv_M * v_t;

    price_std = sqrt(exp(-2 * r * T) * (inv_M * price_std_dev - pow(inv_M * v_t, 2))) / sqrt(M);

    pnl_mat_free(&matrix);
}

void MonteCarlo::price(double t, double &price, double &price_std, const PnlMat *Past, PnlMat *matrix)
{

    int D = this->option->option_size;
    double r = this->model->interest_rate;
    int M = this->sample_number;
    int N = this->fixing_dates_number;
    double T = this->option->maturity;

    double v_t = 0.0;
    double price_std_dev = 0.0;

    pnl_mat_resize(matrix, N + 1, D);

    for (int i = 0; i < M; i++)
    {
        // get_matrix_of_sim(t , matrix);
        this->model->asset(Past, t, matrix, this->rng);
        double phi_j = this->option->payOff(matrix);
        v_t += phi_j;
        price_std_dev += pow(phi_j, 2);
    }

    double inv_M = 1.0 / (double)M;

    price = std::exp(-r * (T - t)) * inv_M * v_t;

    price_std = sqrt(exp(-2 * r * T) * (inv_M * price_std_dev - pow(inv_M * v_t, 2))) / sqrt(M);

    // pnl_mat_free(&matrix);
}

void MonteCarlo::MonteCarlo::get_delta(PnlVect *Vect)
{
    return;
}

void MonteCarlo::calculPAndL(PnlMat *market_data, double &p_and_l)
{
    double step = option->maturity / (double)fixing_dates_number; // step = T/H
    double r = model->interest_rate;
    double p;
    double price_stdev;
    double v;


    PnlMat *past = pnl_mat_new();
    PnlVect *St_i = pnl_vect_new();
    PnlVect *delta = pnl_vect_create(model->model_size);   // deltai
    PnlVect *delta_1 = pnl_vect_create(model->model_size); // delta{i+1}
    PnlVect *calcul = pnl_vect_create(model->model_size);  // pour contenir le la difference
    PnlMat *matrix = pnl_mat_new();
    PnlVect *delta_stdev = pnl_vect_new();
    //
    // traitement de ti = 0
    // price(delta, delta_stdev)
    get_delta(delta);
    price(p, price_stdev);

    v = p - pnl_vect_scalar_prod(delta, model->spots);

    for (int i = 1; i < hedging_date_number + 1; i++)
    {
        double t =i*step;
        pnl_mat_get_row(St_i, market_data, i); // sti
        // calcul de pi<<<<<<<<
        get_cotations(t, past, market_data);
        // // b contient pi
        // price(t, p, price_stdev, Past, matrix);


        // delta_i
        // delta(past ,delta , delta_stdev , t);
        get_delta(delta_1);
        pnl_vect_minus_vect(delta, delta_1); // delta{i-1} - delta{i}
        v = v * exp(r * step) + pnl_vect_scalar_prod(delta, St_i);
        pnl_vect_clone(delta, delta_1);
    }

    // l'erreur de couverture
    p_and_l = v + pnl_vect_scalar_prod(delta, St_i) - option->payOff(market_data);

    // free
    pnl_mat_free(&past);
    pnl_vect_free(&St_i);
    pnl_vect_free(&delta);
    pnl_vect_free(&delta_1);
    pnl_vect_free(&calcul);
    pnl_vect_free(&delta_stdev);
    pnl_mat_free(&matrix);
}

// void MonteCarlo::delta(PnlVect* deltas_vect, PnlVect* stddev_deltas_vect)
// {
//     int D = this->option->option_size;
//     int M = this->sample_number;
//     int N = this->fixing_dates_number;
//     double h = this->fd_step;
//     double r = this->model->interest_rate;
//     double T = this->option->maturity;

//     PnlMat *original_assets = pnl_mat_create(N + 1, D);
//     this->model->asset(original_assets, this->rng);

//     PnlMat *mat_asset_plus = pnl_mat_create(N + 1, D);
//     PnlMat *mat_asset_minus = pnl_mat_create(N + 1, D);

//     PnlVect *st_0 = pnl_vect_create(D);
//     pnl_mat_get_row(st_0, original_assets, 0); // Initial spot prices (S0)

//     // Temporary vectors to store delta values for each simulation
//     PnlVect *deltas_sum = pnl_vect_create(D);
//     PnlVect *deltas_square_sum = pnl_vect_create(D);

//     // Loop over dimensions (assets)
//     for (int d = 0; d < D; d++) {
//         double delta_sum = 0.0;
//         double delta_square_sum = 0.0;

//         // Monte Carlo simulations
//         for (int j = 0; j < M; j++) {
//             this->model->shift_asset(mat_asset_plus, d, h, original_assets);
//             this->model->shift_asset(mat_asset_minus, d, -h, original_assets);

//             double payoff_plus = this->option->payOff(mat_asset_plus);
//             double payoff_minus = this->option->payOff(mat_asset_minus);

//             double delta_j = (payoff_plus - payoff_minus) / (2 * h);

//             // Accumulate delta and squared delta for variance calculation
//             delta_sum += delta_j;
//             delta_square_sum += delta_j * delta_j;
//         }

//         // Calculate the mean delta
//         double delta_mean = delta_sum * exp(-r * T) / (M * pnl_vect_get(st_0, d));
//         pnl_vect_set(deltas_vect, d, delta_mean);

//         // Calculate the standard deviation of delta
//         double delta_var = (delta_square_sum / M - pow(delta_sum / M, 2)) * exp(-2 * r * T) / pow(pnl_vect_get(st_0, d), 2);
//         double delta_stddev = sqrt(delta_var);
//         pnl_vect_set(stddev_deltas_vect, d, delta_stddev);
//     }

//     // Free allocated memory
//     pnl_mat_free(&original_assets);
//     pnl_mat_free(&mat_asset_plus);
//     pnl_mat_free(&mat_asset_minus);
//     pnl_vect_free(&st_0);
//     pnl_vect_free(&deltas_sum);
//     pnl_vect_free(&deltas_square_sum);
// }


// PnlVect* MonteCarlo::delta(const PnlMat *past, PnlVect* deltas_vect, PnlVect* stddev_deltas_vect, double t)
// {
//     /*
//     Usage of the function
//     BEFORE CALLING THE FUNCTION
//     PnlVect* deltas_vect = pnl_vect_new();
//     PnlMat *past = pnl_mat_create(N + 1, D);
//     PnlVect* stddev_deltas_vect = pnl_vect_new();

//     delta(deltas_vect);

//     AFTER CALLING THE FUNCTION
//     pnl_vect_free(&deltas_vect);
//     pnl_vect_free(&stddev_deltas_vect);
//     pnl_mat_free(&past);
//     */
//     int D = this->option->option_size;   // Nb assets
//     int M = this->sample_number;         // Nb Monte Carlo simulations
//     int N = this->fixing_dates_number;   // Nb fixing dates
//     double h = this->fd_step;            // Finite difference step size
//     double r = this->model->interest_rate;
//     double T = this->option->maturity;
    
//     // last index based on time t
//     int index = compute_last_index(t, T, N);

//     // Create original asset matrix
//     PnlMat *original_assets = pnl_mat_create(N + 1, D);
//     this->model->asset(past, t, original_assets, this->rng);

//     // Create matrices for shifted assets
//     PnlMat *mat_asset_plus = pnl_mat_create(N + 1, D);
//     PnlMat *mat_asset_minus = pnl_mat_create(N + 1, D);
//     PnlVect *st_i = pnl_vect_create(D);  // Holds the asset values at time t
//     pnl_mat_get_row(st_i, original_assets, index);

//     // Temporary vectors to store sums and squared sums for each delta
//     PnlVect *deltas_sum = pnl_vect_create(D);
//     PnlVect *deltas_square_sum = pnl_vect_create(D);

//     // Loop over dimensions (assets)
//     for (int d = 0; d < D; d++)
//     {
//         double delta_sum = 0.0;
//         double delta_square_sum = 0.0;

//         // Monte Carlo simulations
//         for (int j = 0; j < M; j++)
//         {
//             // Shift the asset for finite difference calculation
//             this->model->shift_asset(mat_asset_plus, d, t, h, original_assets);
//             this->model->shift_asset(mat_asset_minus, d, t, -h, original_assets);

//             // Calculate payoffs
//             double payoff_plus = this->option->payOff(mat_asset_plus);
//             double payoff_minus = this->option->payOff(mat_asset_minus);

//             // Calculate the delta for the current simulation
//             double delta_j = (payoff_plus - payoff_minus) / (2 * h);

//             // Accumulate the sum of deltas and squared deltas
//             delta_sum += delta_j;
//             delta_square_sum += delta_j * delta_j;
//         }

//         // mean delta
//         double delta_mean = delta_sum * exp(-r * (T - t)) / (M * pnl_vect_get(st_i, d));
//         pnl_vect_set(deltas_vect, d, delta_mean);

//         // variance and standard deviation
//         double delta_var = (delta_square_sum / M - pow(delta_sum / M, 2)) * exp(-2 * r * (T - t)) / pow(pnl_vect_get(st_i, d), 2);
//         double delta_stddev = sqrt(delta_var);
//         pnl_vect_set(stddev_deltas_vect, d, delta_stddev);
//     }

//     // Free memory
//     pnl_mat_free(&original_assets);
//     pnl_mat_free(&mat_asset_plus);
//     pnl_mat_free(&mat_asset_minus);
//     pnl_vect_free(&st_i);
//     pnl_vect_free(&deltas_sum);
//     pnl_vect_free(&deltas_square_sum);
// }


// void MonteCarlo::get_matrix_of_sim(double t , PnlMat *matrix)
// {

//     double T = this->option->maturity;
//     int N = this->fixing_dates_number;
//     int index = compute_last_index(t, T, N);
//     int D = this->option->option_size;

//     PnlVect *dates = pnl_vect_new();
//     get_all_dates(dates, t, index);

//     PnlMat *matrix_sim = pnl_mat_create(D, N - index);

//     if (t == 0.0)
//     {
//         this->model->asset(t, dates, matrix_sim);
//         pnl_mat_set_col(matrix , this->model->spots , 0);
//         pnl_mat_set_subblock(matrix, matrix_sim, 0, index + 1);

//     }
//     else
//     {
//         PnlVect *s_t = pnl_vect_new();
//         PnlMat *cots = pnl_mat_new();
//         get_cotations(t, cots, s_t);

//         PnlVect *col = pnl_vect_create(D);

//         this->model->asset(t, dates, matrix_sim);

//         for (int j = 0; j < N - index; j++)
//         {
//             pnl_mat_get_col(col, matrix_sim, j);
//             pnl_vect_mult_vect_term(col, s_t);
//             pnl_mat_set_col(matrix_sim, col, j);
//         }

//         pnl_mat_set_subblock(matrix, cots, 0, 0);
//         pnl_mat_set_subblock(matrix, matrix_sim, 0, index + 1);

//         pnl_vect_free(&col);
//         pnl_vect_free(&s_t);
//         pnl_mat_free(&cots);
//     }

//     pnl_mat_free(&matrix_sim);
//     pnl_vect_free(&dates);

// }

void MonteCarlo::get_cotations(double t, PnlMat *cots, PnlMat *market_data)
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

    int H = hedging_date_number;
    int N = fixing_dates_number;
    double T = option->maturity;
    int D = option->option_size;

    int i = compute_last_index(t, T, N);

    pnl_mat_resize(cots, i + 2, D);

    PnlVect *col = pnl_vect_create(D);

    for (int j = 0; j < i + 1; j++)
    {
        pnl_mat_get_row(col, market_data, j * H / N);
        pnl_mat_set_row(cots, col, j);
    }

    // s_t
    pnl_mat_get_row(col, market_data, t * H / T);
    pnl_mat_set_row(cots, col, i + 1);

    // free
    pnl_vect_free(&col);
}
