#include "MathUtilsTest.h"
#include "math/MathUtils.h"
#include <sstream>

using namespace rql;

TEST_F(MathUtilsTest, Round)
{
	EXPECT_EQ(3.0, math::round(3.1));
	EXPECT_EQ(4.0, math::round(3.9));
	EXPECT_EQ(4.0, math::round(4.1));
	EXPECT_EQ(5.0, math::round(4.99));
	EXPECT_EQ(4.0, math::round(4.0));

	EXPECT_EQ(-3.0, math::round(-3.1));
	EXPECT_EQ(-4.0, math::round(-3.9));
	EXPECT_EQ(-4.0, math::round(-4.1));
	EXPECT_EQ(-5.0, math::round(-4.99));
	EXPECT_EQ(-4.0, math::round(-4.0));

	EXPECT_EQ(1E10, math::round(1E10+0.1));
}

TEST_F(MathUtilsTest, Trunc)
{
	EXPECT_EQ(3.0, math::trunc(3.1));
	EXPECT_EQ(3.0, math::trunc(3.9));
	EXPECT_EQ(4.0, math::trunc(4.1));
	EXPECT_EQ(4.0, math::trunc(4.99));
	EXPECT_EQ(4.0, math::trunc(4.0));

	EXPECT_EQ(-3.0, math::trunc(-3.1));
	EXPECT_EQ(-3.0, math::trunc(-3.9));
	EXPECT_EQ(-4.0, math::trunc(-4.1));
	EXPECT_EQ(-4.0, math::trunc(-4.99));
	EXPECT_EQ(-4.0, math::trunc(-4.0));

	EXPECT_EQ(1E10, math::trunc(1E10+0.1));
}

TEST_F(MathUtilsTest, Stringify)
{
	EXPECT_EQ("-1", math::stringify<int>(-1));
	EXPECT_EQ("-1", math::stringify<double>(-1));
	EXPECT_EQ("0.1", math::stringify<double>(0.1));
	EXPECT_EQ("0.1", math::stringify(0.1));
}

TEST_F(MathUtilsTest, Finite)
{
	EXPECT_TRUE(math::finite<double>(100.43));
	EXPECT_FALSE(math::finite<double>(std::numeric_limits<double>::infinity()));
	EXPECT_FALSE(math::finite<double>(-std::numeric_limits<double>::infinity()));
	EXPECT_FALSE(math::finite<double>(std::numeric_limits<double>::quiet_NaN()));
	EXPECT_FALSE(math::finite<double>(std::numeric_limits<double>::signaling_NaN()));
	EXPECT_TRUE(math::finite<int>(104));
}

TEST_F(MathUtilsTest, LoadMatrix)
{
	std::stringstream ss;
	ss << "1.0 -1.0\n";
	ss << "-0.5\t0.25\n";
	Eigen::MatrixXd m;
	math::load_matrix_from_stream<double>(ss, m);
	ASSERT_EQ(2, m.rows());
	ASSERT_EQ(2, m.cols());
	ASSERT_EQ(1.0, m(0, 0));
	ASSERT_EQ(-1.0, m(0, 1));
	ASSERT_EQ(-0.5, m(1, 0));
	ASSERT_EQ(0.25, m(1, 1));
}


