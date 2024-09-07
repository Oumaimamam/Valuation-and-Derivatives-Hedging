#ifndef CONVERT_HPP
#define CONVERT_HPP

#include <nlohmann/json.hpp>
#include "../Option/option.hpp"
#include "../black_scholes_model/black_scholes_model.hpp"
#include "../monte_carlo/monte_carlo.hpp";

void convert_json_to_option(Option *ptr_option, nlohmann::json json);

void convert_json_to_model(BlackScholesModel *model, nlohmann::json json);

void convert_json_to_monte_carlo(MonteCarlo *monteCarlo, nlohmann::json json);

#endif