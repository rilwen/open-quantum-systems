#include <gtest/gtest.h>
#include <limits>
#include <vector>
#include "math/SegmentSearch.h"

using namespace rql;

class SegmentSearchTest: public testing::Test
{
};

TEST_F(SegmentSearchTest,LeftInclusive)
{
	std::vector<double> v1(4);
	v1[0] = -0.3;
	v1[1] = 1.4;
	v1[2] = 4.2;
	v1[3] = std::numeric_limits<double>::infinity();
	double v2[4] = {-std::numeric_limits<double>::infinity(), 0.5, 1.5, 2.5};

	EXPECT_EQ(-1, SegmentSearch::binary_search_left_inclusive(v1, v1.size(), -2));
	EXPECT_EQ(-1, SegmentSearch::binary_search_left_inclusive(v1, v1.size(), -std::numeric_limits<double>::infinity()));
	EXPECT_EQ(-1, SegmentSearch::binary_search_left_inclusive(v1, -2));
	EXPECT_EQ(1, SegmentSearch::binary_search_left_inclusive(v1, 2.1));
	for (unsigned int i = 0; i < v1.size(); ++i)
	{
		EXPECT_EQ(i, SegmentSearch::binary_search_left_inclusive(v1, v1[i]));
		EXPECT_EQ(i, SegmentSearch::binary_search_left_inclusive(v1, v1.size(), v1[i]));
	}	

	for (unsigned int i = 0; i < 4; ++i)
	{
		EXPECT_EQ(i, SegmentSearch::binary_search_left_inclusive(v2, 4, v2[i]));
		EXPECT_EQ(i, SegmentSearch::binary_search_left_inclusive(&v2[0], 4, v2[i]));
	}
	EXPECT_EQ(3, SegmentSearch::binary_search_left_inclusive(v2, 4, std::numeric_limits<double>::infinity()));
}

TEST_F(SegmentSearchTest,RightInclusive)
{
	std::vector<double> v1(4);
	v1[0] = -0.3;
	v1[1] = 1.4;
	v1[2] = 4.2;
	v1[3] = std::numeric_limits<double>::infinity();
	double v2[4] = {-std::numeric_limits<double>::infinity(), 0.5, 1.5, 2.5};

	EXPECT_EQ(-1, SegmentSearch::binary_search_right_inclusive(v1, v1.size(), -2));
	EXPECT_EQ(-1, SegmentSearch::binary_search_right_inclusive(v1, -2));
	EXPECT_EQ(1, SegmentSearch::binary_search_right_inclusive(v1, 2.1));
	for (unsigned int i = 0; i < v1.size(); ++i)
	{
		EXPECT_EQ(i-1, SegmentSearch::binary_search_right_inclusive(v1, v1[i]));
		EXPECT_EQ(i-1, SegmentSearch::binary_search_right_inclusive(v1, v1.size(), v1[i]));
	}	

	for (unsigned int i = 0; i < 4; ++i)
	{
		EXPECT_EQ(i-1, SegmentSearch::binary_search_right_inclusive(v2, 4, v2[i]));
		EXPECT_EQ(i-1, SegmentSearch::binary_search_right_inclusive(&v2[0], 4, v2[i]));
	}
}

