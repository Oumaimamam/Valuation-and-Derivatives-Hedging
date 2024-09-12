#include "black_scholes_model.hpp"
#include <cmath>
#include <random>
#include "../Utils/utils.hpp"

BlackScholesModel::BlackScholesModel()
{
    this->volatility = pnl_vect_new();
    this->spots = pnl_vect_new();
}

BlackScholesModel::BlackScholesModel(double rate, PnlVect *vol, PnlVect *spots, double corr, double H, double time_step)
    : interest_rate(rate),
      volatility(vol),
      spots(spots),
      correlation(corr),
      hedging_dates_number(H),
      time_step(time_step)
{
    this->model_size = spots->size;
}

BlackScholesModel::~BlackScholesModel()
{
    pnl_vect_free(&volatility);
    pnl_vect_free(&spots);
}

void BlackScholesModel::asset(const PnlMat *past, double t, PnlMat *path, PnlRng *rng)
{
    int n = path->m; // n = N+1
    int D = this->model_size;
    // matrice L :
    PnlMat *L = pnl_mat_new();
    get_matrix_Cholesky_corr(L);

    PnlVect *L_d = pnl_vect_create(D);
    PnlVect *G = pnl_vect_create(D);
    // t>0
    // col contient le vecteur s_t
    PnlVect *col = pnl_vect_create(D);
    pnl_mat_get_row(col, past, past->m - 1);

    int last_index_t = past->m - 2;

    PnlMat *sub_block = pnl_mat_new();

    pnl_mat_extract_subblock(sub_block, past, 0, last_index_t + 1, 0, D);
    pnl_mat_set_subblock(path, sub_block, 0, 0);
    // calcul de St_i+1
    double r = this->interest_rate;

    PnlVect *calcul = pnl_vect_create(D);

    for (int ligne = last_index_t; ligne < n; ligne++)
    {

        pnl_vect_rng_normal(G, D, rng);
        for (int d = 0; d < D; d++)
        {
            double s_t_i = pnl_vect_get(col, d);
            double sigma_d = pnl_vect_get(this->volatility, d);
            pnl_mat_get_row(L_d, L, d);

            double x = s_t_i * exp((r - pow(sigma_d, 2) / 2.0) * (time_step) + sigma_d * sqrt(time_step) * pnl_vect_scalar_prod(L_d, G)); // t{i+1} - t{i} = T/N
            pnl_vect_set(calcul, d, x);
        }

        pnl_mat_set_row(path, calcul, ligne);
    }

    pnl_vect_free(&calcul);

    // free
    pnl_vect_free(&col);
    pnl_mat_free(&L);
    pnl_vect_free(&L_d);
    pnl_vect_free(&G);
}

// pour le cas t = 0
void BlackScholesModel::asset(PnlMat *path, PnlRng *rng)
{
    int n = path->m; // n = N+1
    int D = this->model_size;
    // matrice L :
    PnlMat *L = pnl_mat_new();
    get_matrix_Cholesky_corr(L);

    PnlVect *L_d = pnl_vect_create(D);
    PnlVect *G = pnl_vect_create(D);
    // t>0
    // col contient le vecteur s_0
    PnlVect *col = pnl_vect_create(D);
    pnl_vect_clone(col, spots);

    pnl_mat_set_row(path, spots, 0);

    // calcul de St_i+1

    double r = this->interest_rate;

    for (int ligne = 1; ligne < n; ligne++)
    {

        pnl_vect_rng_normal(G, D, rng);

        for (int d = 0; d < D; d++)
        {
            double s_t_i = pnl_vect_get(col, d);
            double sigma_d = pnl_vect_get(this->volatility, d);
            pnl_mat_get_row(L_d, L, d);

            double x = s_t_i * exp((r - pow(sigma_d, 2) / 2.0) * (time_step) + sigma_d * sqrt(time_step) * pnl_vect_scalar_prod(L_d, G)); // t{i+1} - t{i} = T/N
            pnl_vect_set(col, d, x);
        }

        pnl_mat_set_row(path, col, ligne);
    }

    // free
    pnl_vect_free(&col);
    pnl_mat_free(&L);
    pnl_vect_free(&L_d);
    pnl_vect_free(&G);
}

void BlackScholesModel::get_matrix_Cholesky_corr(PnlMat *matrix_chol)
{
    /*
    To use this function there are some steps to do :
    PnlMat* matrix = pnl_mat_new();
    get_matrix_Cholesky_corr(matrix);
    pnl_mat_free(&matrix);

    */
    pnl_mat_resize(matrix_chol, this->model_size, this->model_size);
    pnl_mat_set_all(matrix_chol, this->correlation);

    for (int i = 0; i < this->model_size; i++)
        pnl_mat_set_diag(matrix_chol, 1.0, i);

    pnl_mat_chol(matrix_chol);
}

// void BlackScholesModel::shift_asset(PnlMat *shifted_paths, double h, const PnlMat *original_paths)
// {
//     /*
//     To use this function there are some steps to do :
//     PnlMat* shifted_paths = pnl_mat_new();
//     shift_asset(matrix);
//     pnl_mat_free(&shifted_paths);

//     */
//     int D = this->model_size;
//     int N = original_paths->n;

//     pnl_mat_clone(shifted_paths, original_paths);
//     pnl_vect_set( pnl_get_row(original_paths, 0), d, (1 + h ) * pnl_vect_get(pnl_get_row(shifted_paths, 0), d));
//     PnlVect *st_0 =  pnl_get_row(original_paths, 0);

//     for(int i = 1 ; i < N; i++)
//     {
//         pnl_vect_mult_vect_term(pnl_get_row(shifted_paths, i),st_0);

//     }
//     pnl_vect_free(&st_0);

// }

// void BlackScholesModel::shift_asset(PnlMat *shifted_paths, int d,double t, double h, const PnlMat *original_paths)
// {
//     /*
//     To use this function there are some steps to do :
//     PnlMat* shifted_paths = pnl_mat_new();
//     shift_asset(matrix);
//     pnl_mat_free(&shifted_paths);
//     */

//     int D = this->model_size;
//     int N = original_paths->n;
//     int T = this->time_step * N;
//     int index = compute_last_index(t, T, N);

//     pnl_mat_clone(shifted_paths, original_paths);
//     pnl_vect_set( pnl_get_row(original_paths, index), d, (1 + h ) * pnl_vect_get(pnl_get_row(shifted_paths, index), d));
//     PnlVect *st_i =  pnl_get_row(original_paths, index);

//     for(int i = index + 1 ; i < N; i++)
//     {
//         pnl_vect_mult_vect_term(pnl_get_row(shifted_paths, i),st_i);

//     }
//     pnl_vect_free(&st_i);
// }
