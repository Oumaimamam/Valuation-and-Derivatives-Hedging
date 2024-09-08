#include "black_scholes_model.hpp"
#include <cmath>
#include <random>

BlackScholesModel::BlackScholesModel()
{
    this->volatility = pnl_vect_new();
    this->spots = pnl_vect_new();
}

BlackScholesModel::BlackScholesModel(double rate, PnlVect *vol, PnlVect *spots, double corr)
    : interest_rate(rate),
      volatility(vol),
      spots(spots),
      correlation(corr)
{
    this->model_size = spots->size;
}

BlackScholesModel::~BlackScholesModel()
{
    pnl_vect_free(&this->volatility);
    pnl_vect_free(&this->spots);
    pnl_rng_free(&rng);
}

void BlackScholesModel::asset(PnlVect *Dates, PnlMat *mat_asset)
{
    // n = N+1
    int n = Dates->size;
    int D = spots->size;
    // mat_asset = pnl_mat_create(D, n);
    // PnlVect *spots = St0
    PnlVect *col = spots;

    // remplir la prémière colone de la matrice par St0
    pnl_mat_set_col(mat_asset,spots,0);

    // Initialiser le générateur de nombres aléatoires
    rng = pnl_rng_create(PNL_RNG_MERSENNE);
    pnl_rng_sseed(rng, time(NULL));

    // remplir la matrice mat_asset
    for (int j = 1; j < n; j++)
    {
        // calcul de Stj,d
        // col = St_i
        for (int d = 0; d < D; d++)
        {

            double s_t_i = pnl_vect_get(col, d);
            double r = this->interest_rate;
            double sigma_d = pnl_vect_get(this->volatility, d);
            double t_j = pnl_vect_get(Dates, j);
            double t_j_1 = pnl_vect_get(Dates, j - 1);

            // vu u'on travaille avec D= 1, alors on a :
            double L_d = 1;
            // simulier une varibale aléatoire centré reduite dans le cas D=1

            // Générer une variable aléatoire centrée réduite
            double G_i = pnl_rng_normal(rng); // Appeler la fonction pour générer une valeur

            double x = s_t_i * exp((r - pow(sigma_d, 2)/2) * (t_j - t_j_1) + sigma_d * sqrt(t_j - t_j_1) * L_d * G_i);
            pnl_vect_set(col, d, x);
            // col =St_i_1
        }
        pnl_mat_set_col(mat_asset,col,j);
        //free
        pnl_vect_free(&col);
    }
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
