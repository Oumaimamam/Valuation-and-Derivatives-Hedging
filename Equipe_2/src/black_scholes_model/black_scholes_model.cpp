#include "black_scholes_model.hpp"
#include "pnl/pnl_matrix.h"

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
}

void BlackScholesModel::asset()
{
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
