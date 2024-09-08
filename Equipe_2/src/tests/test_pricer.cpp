#include <gtest/gtest.h>

TEST(ExampleTest, HandlesTrueAssertions)
{
    EXPECT_TRUE(true);
}

TEST(ExampleTest, HandlesFalseAssertions)
{
    EXPECT_FALSE(false);
}
