#include <gtest/gtest.h>
#include <iostream>
#include "utils.h"

class UtilsTest: public testing::Test
{
};

TEST_F(UtilsTest,StrConvert)
{
	const std::string a("0.145");
	double da = strConvert<double>(a);
	EXPECT_NEAR(0.145, da, 1E-15);
	da = strConvert<double>(a.c_str());
	EXPECT_NEAR(0.145, da, 1E-15);
	const std::string sa = strConvert<std::string>(da);
	EXPECT_EQ(a, sa);
}
