#include "black_scholes_model.hpp"
#include <cmath>
#include <random>
#include "../Utils/utils.hpp"

BlackScholesModel::BlackScholesModel()
{
    this->volatility = pnl_vect_new();
    this->spots = pnl_vect_new();
}

BlackScholesModel::BlackScholesModel(double rate, PnlVect *vol, PnlVect *spots, double corr, double H)
    : interest_rate(rate),
      volatility(vol),
      spots(spots),
      correlation(corr),
      hedging_dates_number(H)
{
    this->model_size = spots->size;
}

BlackScholesModel::~BlackScholesModel()
{
    pnl_vect_free(&this->volatility);
    pnl_vect_free(&this->spots);
}

// void BlackScholesModel::asset(double t, PnlVect *dates, PnlMat *mat_asset)
// {
//     int n = mat_asset->n;
//     int D = this->model_size;
//     PnlVect *col = pnl_vect_create(D);
//     if (t == 0.0)
//     {
//         pnl_vect_clone(col, this->spots);
//     }
//     else
//     {
//         pnl_vect_set_all(col, 1.0);
//     }

//     // matrice L :
//     PnlMat *L = pnl_mat_new();
//     get_matrix_Cholesky_corr(L);

//     PnlVect *L_d = pnl_vect_create(D);
//     PnlVect *G = pnl_vect_create(D);

//     // remplir la matrice mat_asset
//     for (int j = 0; j < n; j++)
//     {
//         for (int d = 0; d < D; d++)
//         {

//             double s_t_i = pnl_vect_get(col, d);
//             double r = this->interest_rate;
//             double sigma_d = pnl_vect_get(this->volatility, d);
//             double t_j = pnl_vect_get(dates, j);
//             double t_j_1 = pnl_vect_get(dates, j + 1);

//             pnl_mat_get_row(L_d, L, d);
//             pnl_vect_rng_normal(G, D, rng);

//             double x = s_t_i * exp((r - pow(sigma_d, 2)/2.0) * (t_j_1 - t_j) + sigma_d * sqrt(t_j_1 - t_j) * pnl_vect_scalar_prod(L_d, G));
//             pnl_vect_set(col, d, x);
//         }
//         pnl_mat_set_col(mat_asset, col, j);
//     }

//     pnl_vect_free(&col);
//     pnl_mat_free(&L);
//     pnl_vect_free(&L_d);
//     pnl_vect_free(&G);
// }

void BlackScholesModel::asset(const PnlMat *past, double t, PnlMat *path, PnlRng *rng, double maturity)
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
    pnl_mat_get_row(col, past,past->m-1);

    int last_index_t = past->m - 2;

    pnl_mat_extract_subblock(path, past,0, last_index_t +1, 0,D);

    // calcul de St_i+1
    if(t==0.0)
    {
        for(int ligne =  last_index_t; ligne<n; ligne++)
        {

            for (int d = 0; d < D; d++)
            {
                double s_t_i = pnl_vect_get(col, d);
                double r = this->interest_rate;
                double sigma_d = pnl_vect_get(this->volatility, d);
                pnl_mat_get_row(L_d, L, d);
                pnl_vect_rng_normal(G, D, rng);

                double x = s_t_i * exp((r - pow(sigma_d, 2)/2.0) * (maturity/(double)(n-1)) + sigma_d * sqrt(maturity/(double)(n-1)) * pnl_vect_scalar_prod(L_d, G));//t{i+1} - t{i} = T/N
                pnl_vect_set(col, d, x);
            }

            pnl_mat_set_row(path ,col,ligne);
        }

    }else{
        PnlVect *calcul = pnl_vect_create(D);

        for(int ligne =  last_index_t; ligne<n; ligne++)
        {

            for (int d = 0; d < D; d++)
            {
                double s_t_i = pnl_vect_get(col, d);
                double r = this->interest_rate;
                double sigma_d = pnl_vect_get(this->volatility, d);
                pnl_mat_get_row(L_d, L, d);
                pnl_vect_rng_normal(G, D, rng);

                double x = s_t_i * exp((r - pow(sigma_d, 2)/2.0) * (maturity/(double)(n-1)) + sigma_d * sqrt(maturity/(double)(n-1)) * pnl_vect_scalar_prod(L_d, G));//t{i+1} - t{i} = T/N
                pnl_vect_set(calcul, d, x);
            }

            pnl_mat_set_row(path ,calcul,ligne);
        }

        pnl_vect_free(&calcul); 
    }

    //free
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

//past =(st0,...sti,st)  ,  path = vect(st0,...sti,st*S(),...): matrice contient la simulation des sous-jacents
// rng : le générateur
// void BlackScholesModel::get_all_dates(PnlVect *vect, double t, int i, int N, int T) const
// {
//     /*
//   To use this function there are some steps to do :
//   PnlVect* vect = pnl_vect_new();
//   get_all_dates(vect);
//   pnl_vect_free(&vect);
//   */

//     int size = N + 1 - i;
//     pnl_vect_resize(vect, size);
//     double T = option->maturity;
//     double N = (double)N;
//     pnl_vect_set(vect, 0, 0);
//     for (int k = 1; k < size; k++)
//     {
//         double u_k = ((double)(k + i) * T )*(1.0/N )- t;
//         pnl_vect_set(vect, k, u_k);
//     }
// }