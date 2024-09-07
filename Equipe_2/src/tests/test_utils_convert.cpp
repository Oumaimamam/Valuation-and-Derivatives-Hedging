#include <iostream>
#include <fstream>
#include "../Utils/convert.cpp"
#include "../Option/option.hpp"
#include "../monte_carlo/monte_carlo.hpp"

using namespace std;

int main(int argc, char *argv[])
{
    if (argc != 2)
    {
        cout << "Taille attendu 2" << endl;
        return 1;
    }

    test_1_asian_json();

    return 0;
}

void test_1_asian_json()
{
    fstream file("../../../data/asian/asian.json");
    nlohmann::json json = nlohmann::json::parse(file);
    MonteCarlo *monte_carlo = convert_json_to_monte_carlo(json);

    assert(monte_carlo->option->option_size == 2);
    assert(monte_carlo->option->strike == 100.0);

    double spots_expected[2] = {100.0, 100.0};
    assert(is_equals_array(monte_carlo->model->spots->array, spots_expected, 2) == true);

    assert(monte_carlo->option->maturity == 1.5);

    double vols_expected[2] = {0.2, 0.2};

    assert(is_equals_array(monte_carlo->model->volatility->array, vols_expected, 2) == true);

    assert(monte_carlo->model->interest_rate == 0.02);
    assert(monte_carlo->model->correlation == 0.0);

    assert(monte_carlo->option->type == OptionType::Asian);
    double payoff_coff_expected[2] = {0.5, 0.5};
    assert(is_equals_array(monte_carlo->option->payoff_coeffcients->array, payoff_coff_expected, 2) == true);

    assert(monte_carlo->fixing_dates_number == 24);
    assert(monte_carlo->sample_number == 50000);

    file.close();
}

bool is_equals_array(double *ptr_arr_res, double *arr_expected, int size)
{
    // Comparaison manuelle élément par élément
    for (int i = 0; i < size; ++i)
    {
        if (ptr_arr_res[i] == arr_expected[i])
        {
            return false;
        }
    }
    return true;
}