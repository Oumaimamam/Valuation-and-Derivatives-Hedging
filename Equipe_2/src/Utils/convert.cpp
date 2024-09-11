#include "convert.hpp"
#include "../json_helper.hpp"
#include <iostream>
#include <fstream>
#include "../Option/option_type.hpp"
#include "../Option/option.hpp"
#include "../black_scholes_model/black_scholes_model.hpp"
#include "../monte_carlo/monte_carlo.hpp"

Option *convert_json_to_option(nlohmann::json json)
{
    double T;
    double K;
    PnlVect *coeff;
    double D;
    std::string type_str;

    json.at("maturity").get_to(T);

    if (json.contains("strike"))
    {

        json.at("strike").get_to(K);
    }
    else
    {
        K = 0.0;
    }
    json.at("option size").get_to(D);
    json.at("payoff coefficients").get_to(coeff);
    if (coeff->size == 1 && D > 1)
    {
        pnl_vect_resize_from_scalar(coeff, D, GET(coeff, 0));
    }
    json.at("option type").get_to(type_str);

    OptionType type = stringToOptionType(type_str);

    // Option* option = new Option();
    Option *option = Option::GetOption(T, K, type, D, coeff);

    return option;
}

BlackScholesModel *convert_json_to_model(nlohmann::json json)
{
    double r;
    double size;
    PnlVect *vols;
    PnlVect *spots;
    double corr;
    double H;

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
    json.at("hedging dates number").get_to(H);

    BlackScholesModel *model = new BlackScholesModel(r, vols, spots, corr, H);

    return model;
}

MonteCarlo *convert_json_to_monte_carlo(std::string file_path, std::string data_file_name)
{
    int N;
    int M;
    double h;

    std::ifstream file(file_path);
    if (!file.is_open())
    {
        std::cerr << "Error opening file" << std::endl;
        exit(1);
    }
    nlohmann::json json = nlohmann::json::parse(file);

    json.at("fixing dates number").get_to(N);
    json.at("sample number").get_to(M);
    json.at("fd step").get_to(h);

    BlackScholesModel *model = convert_json_to_model(json);
    Option *option = convert_json_to_option(json);

    file.close();

    if (data_file_name != "")
    {
        PnlMat *data = pnl_mat_create_from_file(data_file_name.c_str());
        MonteCarlo *monte_carlo = new MonteCarlo(option, model, N, M, data, h);

        return monte_carlo;
    }
    else
    {
        MonteCarlo *monte_carlo = new MonteCarlo(option, model, N, M, NULL, h);
        return monte_carlo;
    }
}