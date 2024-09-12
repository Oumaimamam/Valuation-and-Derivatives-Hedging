#include "monte_carlo.hpp"
#include <cmath>
#include "../Utils/convert.hpp"
// #include "../Utils/construct_append_mat.hpp"
#include "../Utils/utils.hpp"

// void MonteCarlo::get_all_dates(PnlVect *vect, double t, int i) const
// {
//     /*
//   To use this function there are some steps to do :
//   PnlVect* vect = pnl_vect_new();
//   get_all_dates(vect);
//   pnl_vect_free(&vect);

//   */
//     int size = this->fixing_dates_number + 1 - i;
//     pnl_vect_resize(vect, size);
//     double T = option->maturity;
//     double N = (double)this->fixing_dates_number;
//     pnl_vect_set(vect, 0, 0);
//     for (int k = 1; k < size; k++)
//     {
//         double u_k = ((double)(k + i) * T )*(1.0/N )- t;
//         pnl_vect_set(vect, k, u_k);
//     }
// }

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
    // Libérer la mémoire si nécessaire
    if (option != nullptr)
    {
        delete option;
    }
    if (model != nullptr)
    {
        delete model;
    }


    if (rng != nullptr)
    {
        pnl_rng_free(&rng);
    }
    

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
        // get_matrix_of_sim(t , matrix);
        this->model->asset(matrix,this->rng);
        double phi_j = this->option->payOff(matrix);
        v_t += phi_j;
        price_std_dev += pow(phi_j, 2);
    }


    double inv_M = 1.0 / (double)M;

    price = std::exp(-r * (T )) * inv_M * v_t;

    price_std = sqrt(exp(-2 * r * T) * (inv_M * price_std_dev - pow(inv_M * v_t, 2))) / sqrt(M);

    pnl_mat_free(&matrix);
}


void MonteCarlo::price(double t, double &price, double &price_std, const PnlMat *Past)
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
        // get_matrix_of_sim(t , matrix);
        this->model->asset(Past,t,matrix,this->rng);
     double phi_j = this->option->payOff(matrix);
        v_t += phi_j;
        price_std_dev += pow(phi_j, 2);
    }


    double inv_M = 1.0 / (double)M;

    price = std::exp(-r * (T - t)) * inv_M * v_t;

    price_std = sqrt(exp(-2 * r * T) * (inv_M * price_std_dev - pow(inv_M * v_t, 2))) / sqrt(M);

    pnl_mat_free(&matrix);
}



PnlVect* MonteCarlo::delta(PnlVect* deltas_vect, double t, const PnlMat *cots, PnlVect *s_t)
{
    // Usage of the function
    // BEFORE CALLING THE FUNCTION
    //PnlVect* deltas_vect = pnl_vect_new();
    //PnlVect* s_t = pnl_vect_new();
    //PnlMat* cots = pnl_mat_new();

    //get_cotations(t, cots, s_t);

    // AFTER CALLING THE FUNCTION
    //delta(deltas_vect, t, cots, s_t);
    //pnl_vect_free(&deltas_vect);
    //pnl_vect_free(&s_t);
    //pnl_mat_free(&cots);


    double T = this->option->maturity;
    int M = this->sample_number; 
    double h = this->fd_step;
    int i = compute_last_index(t, T, N);

    PnlMat *mat_asset = pnl_mat_create(N - i, D);
    PnlMat *mat_asset_plus = pnl_mat_new(); 
    PnlMat *mat_asset_minus = pnl_mat_new();
    PnlVect *gd_plus = pnl_vect_new();
    PnlVect *gd_minus = pnl_vect_new();
    PnlMat *M_plus = pnl_mat_create(N + 1, D);
    PnlMat *M_minus = pnl_mat_create(N + 1, D);


    for (int d = 0; d < D; d++)
    {
        double delta_sum = 0.0;

        for (int j = 0; j < M; j++)
        {

            // gd+ and gd-
            pnl_vect_clone(gd_plus, s_t);
            pnl_vect_set(gd_plus, d, pnl_vect_get(s_t, d) * (1 + h));
            pnl_vect_clone(gd_minus, s_t);
            pnl_vect_set(gd_minus, d, pnl_vect_get(s_t, d) * (1 - h));

            //asset(cots, t, PnlMat *path, this->rng, this->option->maturity)
            //get_all_dates(dates);
            //this->model->asset(t, dates, mat_asset);

            // matrix of sim
            get_matrix_of_sim(t, mat_asset);
            pnl_mat_clone(mat_asset_plus, mat_asset);
            pnl_mat_clone(mat_asset_minus, mat_asset);

            for ( int p = 0; p < N ; p++)
            {
                pnl_vect_mult_vect_term(mat_asset_plus[p], gd_plus[p]);
                pnl_vect_mult_vect_term(mat_asset_minus[p], gd_minus[p]);
            }

            pnl_mat_set_subblock(M_plus, cots, 0, D);
            pnl_mat_set_subblock(M_plus, mat_asset_plus, i, D);
            pnl_mat_set_subblock(M_minus, cots, 0, D);
            pnl_mat_set_subblock(M_minus, mat_asset_minus, i, D);
            
            double payoff_minus = this->option->payOff(M_plus);
            double payoff_plus = this->option->payOff(M_minus);
            
            delta_sum += (payoff_plus - payoff_minus) / (2 * h);
        }

        
        double delta_d = delta_sum * exp(- r * (T - t)) / (M * pnl_vect_get(s_t, d));
        pnl_vect_set(deltas_vect, d, delta_d);
    }

    pnl_vect_free(&gd_plus);
    pnl_vect_free(&gd_minus);
    pnl_mat_free(&mat_asset);
    pnl_mat_free(&mat_asset_plus);
    pnl_mat_free(&mat_asset_minus);
    pnl_mat_free(&M_plus);
    pnl_mat_free(&M_minus);

    return deltas_vect;
}
PnlVect* MonteCarlo::delta(PnlVect* deltas_vect, double t, const PnlMat *cots, PnlVect *s_t)
{
    // Usage of the function
    // BEFORE CALLING THE FUNCTION
    //PnlVect* deltas_vect = pnl_vect_new();
    //PnlVect* s_t = pnl_vect_new();
    //PnlMat* cots = pnl_mat_new();

    //get_cotations(t, cots, s_t);

    // AFTER CALLING THE FUNCTION
    //delta(deltas_vect, t, cots, s_t);
    //pnl_vect_free(&deltas_vect);
    //pnl_vect_free(&s_t);
    //pnl_mat_free(&cots);


    double T = this->option->maturity;
    int M = this->sample_number; 
    double h = this->fd_step;
    int i = compute_last_index(t, T, N);

    PnlMat *mat_asset = pnl_mat_create(N - i, D);
    PnlMat *mat_asset_plus = pnl_mat_new(); 
    PnlMat *mat_asset_minus = pnl_mat_new();
    PnlVect *gd_plus = pnl_vect_new();
    PnlVect *gd_minus = pnl_vect_new();
    PnlMat *M_plus = pnl_mat_create(N + 1, D);
    PnlMat *M_minus = pnl_mat_create(N + 1, D);


    for (int d = 0; d < D; d++)
    {
        double delta_sum = 0.0;

        for (int j = 0; j < M; j++)
        {

            // gd+ and gd-
            pnl_vect_clone(gd_plus, s_t);
            pnl_vect_set(gd_plus, d, pnl_vect_get(s_t, d) * (1 + h));
            pnl_vect_clone(gd_minus, s_t);
            pnl_vect_set(gd_minus, d, pnl_vect_get(s_t, d) * (1 - h));

            //asset(cots, t, PnlMat *path, this->rng, this->option->maturity)
            //get_all_dates(dates);
            //this->model->asset(t, dates, mat_asset);

            // matrix of sim
            get_matrix_of_sim(t, mat_asset);
            pnl_mat_clone(mat_asset_plus, mat_asset);
            pnl_mat_clone(mat_asset_minus, mat_asset);

            for ( int p = 0; p < N ; p++)
            {
                pnl_vect_mult_vect_term(mat_asset_plus[p], gd_plus[p]);
                pnl_vect_mult_vect_term(mat_asset_minus[p], gd_minus[p]);
            }

            pnl_mat_set_subblock(M_plus, cots, 0, D);
            pnl_mat_set_subblock(M_plus, mat_asset_plus, i, D);
            pnl_mat_set_subblock(M_minus, cots, 0, D);
            pnl_mat_set_subblock(M_minus, mat_asset_minus, i, D);
            
            double payoff_minus = this->option->payOff(M_plus);
            double payoff_plus = this->option->payOff(M_minus);
            
            delta_sum += (payoff_plus - payoff_minus) / (2 * h);
        }

        
        double delta_d = delta_sum * exp(- r * (T - t)) / (M * pnl_vect_get(s_t, d));
        pnl_vect_set(deltas_vect, d, delta_d);
    }

    pnl_vect_free(&gd_plus);
    pnl_vect_free(&gd_minus);
    pnl_mat_free(&mat_asset);
    pnl_mat_free(&mat_asset_plus);
    pnl_mat_free(&mat_asset_minus);
    pnl_mat_free(&M_plus);
    pnl_mat_free(&M_minus);

    return deltas_vect;
}
PnlVect* MonteCarlo::delta(PnlVect* deltas_vect)
{
    /*
    Usage of the function
    BEFORE CALLING THE FUNCTION
    PnlVect* deltas_vect = pnl_vect_new();

    delta(deltas_vect);

    AFTER CALLING THE FUNCTION
    pnl_vect_free(&deltas_vect);
    */

    int D = this->option->option_size;
    int M = this->sample_number;
    int N = this->fixing_dates_number;
    double h = this->fd_step;

    PnlMat *original_assets = pnl_mat_create(N + 1, D);
    this->model->asset(original_assets, this->rng);
    PnlMat *mat_asset_plus = pnl_mat_create(N + 1, D); 
    PnlMat *mat_asset_minus = pnl_mat_create(N + 1, D);
    PnlVect *st_0 =  pnl_get_row(original_assets, 0);
    for (int d = 0; d < D; d++)
    {
        double delta_sum = 0.0;
        for (int j = 0; j < M; j++)
        {
            this->model->shift_asset(mat_asset_plus, h,original_assets);
            this->model->shift_asset(mat_asset_minus, -h,original_assets);
            double payoff_minus = this->option->payOff(mat_asset_plus);
            double payoff_plus = this->option->payOff(mat_asset_minus);
            delta_sum += (payoff_plus - payoff_minus) / (2 * h);
        }

        double delta_d = delta_sum * exp(- r * T ) / (M * pnl_vect_get(st_0, d));
        pnl_vect_set(deltas_vect, d, delta_d);
    }
    pnl_mat_free(&original_assets);
    pnl_mat_free(&mat_asset_plus);
    pnl_mat_free(&mat_asset_minus);
    pnl_vect_free(&st_0);

    return deltas_vect;
}


PnlVect* MonteCarlo::delta(const PnlMat *past, PnlVect* deltas_vect, double t)
{
    /*
    Usage of the function
    BEFORE CALLING THE FUNCTION
    PnlVect* deltas_vect = pnl_vect_new();
    PnlVect* s_t = pnl_vect_new();
    delta(deltas_vect);
    AFTER CALLING THE FUNCTION
    pnl_vect_free(&deltas_vect);
    pnl_vect_free(&s_t);
    */
    int D = this->option->option_size;
    int M = this->sample_number;
    int N = this->fixing_dates_number;
    double h = this->fd_step;
    int T = this->model->time_step * N;
    int index = compute_last_index(t, T, N);
    PnlMat *past = pnl_mat_create(index, D);
    PnlMat *original_assets = pnl_mat_create(N + 1, D);
    this->model->asset(past, t, original_assets,  this->rng);

    PnlMat *mat_asset_plus = pnl_mat_create(N + 1, D); 
    PnlMat *mat_asset_minus = pnl_mat_create(N + 1, D);
    PnlVect *st_i =  pnl_get_row(original_assets, index);

    for (int d = index + 1; d < D; d++)
    {
        double delta_sum = 0.0;
        for (int j = 0; j < M; j++)
        {
            this->model->shift_asset(mat_asset_plus, d, h,original_assets);
            this->model->shift_asset(mat_asset_minus, d, -h,original_assets);
            double payoff_minus = this->option->payOff(mat_asset_plus);
            double payoff_plus = this->option->payOff(mat_asset_minus);
            delta_sum += (payoff_plus - payoff_minus) / (2 * h);
        }

        double delta_d = delta_sum * exp(- r * T ) / (M * pnl_vect_get(st_i, d));
        pnl_vect_set(deltas_vect, d, delta_d);
    }
    pnl_mat_free(&original_assets);
    pnl_mat_free(&mat_asset_plus);
    pnl_mat_free(&mat_asset_minus);
    pnl_vect_free(&st_i);
    return deltas_vect; 
}



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











// void MonteCarlo::get_cotations(double t, PnlMat *cots, PnlVect *s_t)
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

//     if (market_data == NULL)
//     {
//         throw std::invalid_argument("argument data.txt non fourni");
//         exit(1);
//     }

//     int H = this->model->hedging_dates_number;
//     int N = this->fixing_dates_number;
//     double T = this->option->maturity;
//     int D = this->option->option_size;

//     int i = compute_last_index(t, T, N);

//     pnl_mat_resize(cots, D, i + 1);

//     PnlVect *col = pnl_vect_create(D);

//     for (int j = 0; j < i + 1; j++)
//     {
//         pnl_mat_get_row(col, this->market_data, j * H / N);
//         pnl_mat_set_col(cots, col, j);
//     }

//     int index_t = t * H / T;
//     pnl_mat_get_row(s_t, this->market_data, index_t);

//     pnl_vect_free(&col);
// }
