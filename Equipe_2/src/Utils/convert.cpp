#include "convert.hpp"
#include "../json_helper.hpp"
Option *convert_json_to_option(nlohmann::json json)
{
    return nullptr;
}

BlackScholesModel *convert_json_to_model(nlohmann::json json)
{
    double r;
    double size;
    PnlVect *vols;
    PnlVect *spots;
    double corr;

    json.at("interest rate").get_to(r);
    json.at("option size").get_to(size);
    json.at("volatility").get_to(vols);
    if (vols->size == 1 && size > 1)
    {
        pnl_vect_resize_from_scalar(vols, size, GET(vols, 0));
    }
    json.at("spot").get_to(spots);
    if (spots->size == 1 && size > 1)
    {
        pnl_vect_resize_from_scalar(spots, size, GET(spots, 0));
    }
    json.at("correlation").get_to(corr);

    BlackScholesModel *model = new BlackScholesModel(r, vols, spots, corr);

    return model;
}

MonteCarlo *convert_json_to_monte_carlo(nlohmann::json json)
{
    int N;
    int M;

    json.at("fixing dates number").get_to(N);
    json.at("sample number").get_to(M);

    BlackScholesModel *model = convert_json_to_model(json);
    Option *option = convert_json_to_option(json);

    MonteCarlo *monte_carlo = new MonteCarlo(option, model, N, M);

    return monte_carlo;
}