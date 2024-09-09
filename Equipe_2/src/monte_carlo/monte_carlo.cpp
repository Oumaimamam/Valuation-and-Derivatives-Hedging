#include "monte_carlo.hpp"
#include <cmath>
#include "../Utils/convert.hpp"
#include "../Utils/utils.hpp"

void MonteCarlo::get_all_dates(PnlVect *vect , double t , int i ) const
{
    /*
  To use this function there are some steps to do :
  PnlVect* vect = pnl_vect_new();
  get_all_dates(vect);
  pnl_vect_free(&vect);

  */
    int size = this->fixing_dates_number + 1 - i;
    pnl_vect_resize(vect, size);
    // Maturity
    double T = option->maturity;

    pnl_vect_set(vect, 0, 0);
    for (int k = 1; k < size; k++)
    {
        double u_k = ((k +i)* T / this->fixing_dates_number) - t;
        pnl_vect_set(vect, k, u_k);
        
    }
}

MonteCarlo::MonteCarlo(Option *option, BlackScholesModel *model, int N, int M)
    : option(option),
      model(model),
      fixing_dates_number(N),
      sample_number(M)
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
}

double MonteCarlo::price(double t)
{
    // calcul du prix à l'instat t = 0 et pour d = 1
    double v_0 = 0.0;

    int D = this->option->option_size;

    int i = compute_last_index(t , this->option->maturity , this->fixing_dates_number);

    PnlVect *dates = pnl_vect_new();
    get_all_dates(dates , t , i);


    PnlVect* s_t = pnl_vect_new();
    PnlMat* cots = pnl_mat_new();


    PnlVect* spots = pnl_vect_create(D);

    if(t == 0.0 ) {
        pnl_vect_clone(spots , this->model->spots);
    } else {
        pnl_vect_set_all(spots , 1.0);

    }


    get_cotations(t , cots , s_t);



    PnlMat *matrix = pnl_mat_create(D, this->fixing_dates_number + 1);

    PnlMat* matrix_sim = pnl_mat_create(D , this->fixing_dates_number - i);

    PnlVect* col = pnl_vect_create(D);


    for (int i = 1; i < this->sample_number + 1; i++)
    {
        this->model->asset(spots , dates, matrix_sim);

        for (int j = 0; j < this->fixing_dates_number - i ; i++)
        {
            pnl_mat_get_col(col,matrix_sim , j);
            pnl_vect_mult_vect_term( col,s_t );  //In place ?? :  col ?? or s_t ??
            pnl_mat_set_col(matrix_sim , col , j);
        }
        
        pnl_mat_set_subblock(matrix , cots , 0 , 0) ;
        pnl_mat_set_subblock(matrix , matrix_sim , 0 , i+1);

        v_0 += this->option->payOff(matrix);

    }

    double r = this->model->interest_rate;
    int T = this->option->maturity;

    // free :
    pnl_vect_free(&dates);
    pnl_vect_free(&col);
    pnl_mat_free(&matrix);
    pnl_vect_free(&s_t);
    pnl_mat_free(&cots);
    pnl_vect_free(&spots);
    pnl_mat_free(&matrix_sim);

    return std::exp(-r * (T - t)) * (1 / this->sample_number) * v_0;
}

void MonteCarlo::get_cotations(double t, PnlMat *cots , PnlVect* s_t)
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

    PnlMat* data = pnl_mat_create_from_file("../../data/call/call_market.txt");
    int H = this->model->hedging_dates_number;
    int N = this->fixing_dates_number;
    double T = this->option->maturity; 
    int D = this->option->option_size;

    int i = compute_last_index(t , T , N);

    pnl_mat_resize(cots , D , i + 1);

    PnlVect* col = pnl_vect_create(D);

    for(int j = 0 ; j < i + 1 ;  j++ ) {
        pnl_mat_get_col(col , data , j*H/N);
        pnl_mat_set_col(cots, col , j );
    }

    int index_t = t*H/T;
    pnl_mat_get_col(s_t , data , index_t);

    pnl_mat_free(&data);
    pnl_vect_free(&col);

}



// Destructeur

// Méthode privée pour obtenir les dates
// PnlVect *MonteCarlo::getDates() const
// {
//     PnlVect *Dates = pnl_vect_create(fixing_dates_number + 1);
//     // pnl_vect_set(PnlVect ∗v, int i, double x) to set x in place i
//     for (int k = 0, k < fixing_dates_number + 1, k++)
//     {
//         pnl_vect_set(Dates, k, k * option.getMaturity() / fixing_dates_number);
//     }
// }

// Méthode pour calculer le prix
// double MonteCarlo::price(double t)
// {
//     double r = model.interesetRate;
//     int T = option.maturity;
//     return exp(-r * T) * (1 / sample_number) * sum();
// }

// double MonteCarlo::sum()
// {
//     double res = 0.0;
//     for (int j = 1, j < fixing_dates_number + 1, j++)
//     {
//         pnlMat *vectSim = pnl_mat_create(1, fixing_dates_number + 1);
//         res += option.payoff(asset(getDates()), vectSim);
//     }
//     return res / sample_number;
// }
