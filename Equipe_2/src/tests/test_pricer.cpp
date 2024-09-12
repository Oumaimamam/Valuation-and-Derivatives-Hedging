#include <gtest/gtest.h>
#include "../Utils/convert.hpp"
#include "../monte_carlo/monte_carlo.hpp"
#include <iostream>

TEST(MonteCarloTest, TestingPriceCall)
{
    MonteCarlo *monte_carlo = convert_json_to_monte_carlo("../../data/call/call.json");
    double price = 0.0;
    double price_std_dev = 0.0;
    monte_carlo->price(0, price, price_std_dev);

    std::cout << "price_t = " << price << std::endl;
    std::cout << "expected_t = " << 10.549999532823513 << std::endl;

    std::cout << "price_std_dev = " << price_std_dev << std::endl;
    std::cout << "expected_price_std_dev = " << 0.06638859132968258 << std::endl;
    delete monte_carlo;
}

TEST(MonteCarloTest, TestingPriceAsian)
{
    MonteCarlo *monte_carlo = convert_json_to_monte_carlo("../../data/asian/asian.json");
    double price = 0.0;
    double price_std_dev = 0.0;
    monte_carlo->price(0, price, price_std_dev);

    std::cout << "price_t = " << price << std::endl;
    std::cout << "expected_t = " << 4.6230359733474184 << std::endl;
    delete monte_carlo;
}

TEST(MonteCarloTest, TestingPriceBasket2D)
{
    MonteCarlo *monte_carlo = convert_json_to_monte_carlo("../../data/basket/basket_2d/basket_2d.json");

    double price = 0.0;
    double price_std_dev = 0.0;

    monte_carlo->price(0, price, price_std_dev);

    std::cout << "price_t = " << price << std::endl;
    std::cout << "expected_t = " << 8.340413215635659 << std::endl;
    delete monte_carlo;
}

TEST(MonteCarloTest, TestingPriceBasket5D)
{
    MonteCarlo *monte_carlo = convert_json_to_monte_carlo("../../data/basket/basket_5d/basket_5d.json");
    double price = 0.0;
    double price_std_dev = 0.0;
    monte_carlo->price(0, price, price_std_dev);
    std::cout << "price_t = " << price << std::endl;
    std::cout << "expected_t = " << 6.362373479357343 << std::endl;
    delete monte_carlo;
}

TEST(MonteCarloTest, TestingPriceBasket5D1)
{
    MonteCarlo *monte_carlo = convert_json_to_monte_carlo("../../data/basket/basket_5d_1/basket_5d_1.json");
    double price = 0.0;
    double price_std_dev = 0.0;
    monte_carlo->price(0, price, price_std_dev);
    std::cout << "price_t = " << price << std::endl;
    std::cout << "expected_t = " << 17.76738571855691 << std::endl;
    delete monte_carlo;
}

TEST(MonteCarloTest, TestingPriceBasket40D)
{
    MonteCarlo *monte_carlo = convert_json_to_monte_carlo("../../data/basket/basket_40d/basket_40d.json");
    double price = 0.0;
    double price_std_dev = 0.0;

    monte_carlo->price(0, price, price_std_dev);
    std::cout << "price_t = " << price << std::endl;
    std::cout << "expected_t = " << 9.1280832137086 << std::endl;
    delete monte_carlo;
}

TEST(MonteCarloTest, TestingPricePerf)
{
    MonteCarlo *monte_carlo = convert_json_to_monte_carlo("../../data/perf/perf.json");
    double price = 0.0;
    double price_std_dev = 0.0;

    monte_carlo->price(0, price, price_std_dev);

    std::cout << "price_t = " << price << std::endl;
    std::cout << "expected_t = " << 1.2575257621859213 << std::endl;
    delete monte_carlo;
}
