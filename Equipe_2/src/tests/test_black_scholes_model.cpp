#include <gtest/gtest.h>
#include "../monte_carlo/monte_carlo.hpp"
#include "../Utils/convert.hpp"

TEST(BlackScholesModel, TestAssetCall)
{
    MonteCarlo* monte_carlo = convert_json_to_monte_carlo("../../data/call/call.json");

    BlackScholesModel* model = monte_carlo->model;

    PnlVect* list_ti = pnl_vect_new();
    monte_carlo->get_all_dates(list_ti);

    PnlMat* mat_asset = pnl_mat_create( monte_carlo->option->option_size,monte_carlo->fixing_dates_number +1 );

    model->asset(list_ti , mat_asset);

    pnl_vect_free(&list_ti);
    pnl_mat_free(&mat_asset);

}


