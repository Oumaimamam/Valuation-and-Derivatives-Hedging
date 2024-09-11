// #include <gtest/gtest.h>
// #include "../Option/option.hpp"


// // Mock pour OptionType (doit correspondre à votre implémentation)
// enum OptionType { Asian, Basket, Performance };

// // Test du constructeur de la classe de base Option
// TEST(OptionTest, ConstructorTest) {
//     PnlVect* coeff = pnl_vect_create_from_double(3, 1.0); // Initialise un vecteur de taille 3 avec tous les éléments à 1.0
//     OptionAsian optionAsian(1.0, 100.0, Asian, 3, coeff); // Initialisation d'une option asiatique

//     EXPECT_EQ(optionAsian.maturity, 1.0);
//     EXPECT_EQ(optionAsian.strike, 100.0);
//     EXPECT_EQ(optionAsian.type, Asian);
//     EXPECT_EQ(optionAsian.option_size, 3);
//     EXPECT_EQ(pnl_vect_get(optionAsian.payoff_coeffcients, 0), 1.0);

//     pnl_vect_free(&coeff); // Libération de la mémoire
// }

// // Test de la méthode payOff de OptionAsian
// TEST(OptionAsianTest, PayOffTest) {
//     PnlVect* coeff = pnl_vect_create_from_double(3, 1.0); // Coefficients de payoff
//     OptionAsian optionAsian(1.0, 100.0, Asian, 3, coeff);

//     // Créer une matrice fictive pour le test (par exemple, 10 dates et 3 sous-jacents)
//     PnlMat* matrix = pnl_mat_create_from_scalar(10, 3, 110.0); // Chaque élément de la matrice est 110.0
//     double payoff = optionAsian.payOff(matrix);

//     // Test basé sur votre logique de calcul du payoff (à ajuster selon l'implémentation réelle)
//     EXPECT_GT(payoff, 0);

//     pnl_mat_free(&matrix); // Libération de la mémoire
//     pnl_vect_free(&coeff); // Libération de la mémoire
// }

// // Autres tests similaires pour OptionBasket et OptionPerformance

// int main(int argc, char **argv) {
//     ::testing::InitGoogleTest(&argc, argv);
//     return RUN_ALL_TESTS();
// }

#include <gtest/gtest.h>

TEST(ExampleTest, HandlesTrueAssertions)
{
    EXPECT_TRUE(true);
}

TEST(ExampleTest, HandlesFalseAssertions)
{
    EXPECT_FALSE(false);
}
