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
}

BlackScholesModel::~BlackScholesModel()
{
    pnl_vect_free(&this->volatility);
    pnl_vect_free(&this->spots);
}

void BlackScholesModel::asset(PnlVect *Dates)
{
    int N = Dates->size+1;
    int D = spots->size;
    mat_asset = pnl_mat_create(1,N);
    //PnlVect *spots = St0
    PnlVect* col =spots ;

    // remplir la prémière colone de la matrice par St0
    for (int i = 0; i < 1; ++i) {
        pnl_mat_set(mat_asset, i, 0, pnl_vect_get(spots, i));  // Copier chaque élément du vecteur dans la première colonne de la matrice
    }


    // remplir la matrice mat_asset
    for(int j= 1 )
    //remplir St_i
    for(int j=1; j<N; j++)
    {
        // calcul de Stj,d
        for(int d = 0; d < D; d++)
        {

            double s_t_i = pnl_vect_get(col,j-1); 
            double r = this->interest_rate;
            double sigma_d = pnl_vect_get(this->volatility , d);
            double t_j = pnl_vect_get(Dates,j);
            double t_j_1 = pnl_vect_get(Dates,j-1);

            // vu u'on travaille avec D= 1, alors on a :
            double L_d = sqrt(this->correlation);
            // simulier une varibale aléatoire centré reduite dans le cas D=1
            std::random_device rd;  // Obtenir une valeur aléatoire non déterministe
            std::mt19937 gen(rd()); // Initialiser le générateur avec la valeur obtenue

            // Créer une distribution normale centrée réduite (moyenne=0, écart-type=1)
            std::normal_distribution<double> d(0.0, 1.0);

            // Générer un échantillon de la distribution
            double G_i = d(gen);



            double x = s_t_i*exp((r - pow(sigma_d,2)) * (t_j - t_j_1) + sigma_d*sqrt(t_j - t_j_1)*L_d* G_i) ;
            pnl_vect_set(col, d, x);
        }

    }
    
}
