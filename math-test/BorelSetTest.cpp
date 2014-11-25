#include "BorelSetTest.h"
#include "math/BorelSet.h"
#include <list>
#include <limits>
#include <iostream>

using namespace rql::math;

TEST_F(BorelSetTest,ConstructionIterators)
{
	std::list<Segment> input;
	input.push_back(Segment(0.1, 0.2, false, false));
	input.push_back(Segment(0.2, 0.25, false, false));
	input.push_back(Segment(0.4, 0.5, false, false));
	const BorelSet bs1(input.begin(), input.end(), true);
	EXPECT_EQ(2, bs1.nbr_segments());
	EXPECT_EQ(Segment(0.1, 0.25, true, true), bs1.segment(0));
	EXPECT_EQ(Segment(0.4, 0.5, true, true), bs1.segment(1));
	const BorelSet bs2(input.begin(), input.end());
	EXPECT_EQ(3, bs2.nbr_segments());
	EXPECT_EQ(Segment(0.1, 0.2, false, false), bs2.segment(0));
	EXPECT_EQ(Segment(0.2, 0.25, false, false), bs2.segment(1));
	EXPECT_EQ(Segment(0.4, 0.5, false, false), bs2.segment(2));
	input.push_back(Segment::EMPTY);
	const BorelSet bs3(input.begin(), input.end());
	EXPECT_EQ(bs2, bs3);
}

TEST_F(BorelSetTest,ConstructionSingle)
{
	const Segment s1(0.1, 0.2, true, false);
	const BorelSet bs1(s1);
	EXPECT_EQ(1, bs1.nbr_segments());
	EXPECT_EQ(s1, bs1.segment(0));
	const BorelSet bs2;
	EXPECT_EQ(1, bs2.nbr_segments());
	EXPECT_EQ(Segment::EMPTY, bs2.segment(0));
}

TEST_F(BorelSetTest,IsEmpty)
{
	const BorelSet bs;
	EXPECT_TRUE(bs.is_empty());
}

TEST_F(BorelSetTest,Comparison)
{
	std::list<Segment> input;
	input.push_back(Segment(0.1, 0.2, false, false));
	input.push_back(Segment(0.2, 0.25, false, false));
	input.push_back(Segment(0.4, 0.5, false, false));
	const BorelSet bs1a(input.begin(), input.end(), true);
	const BorelSet bs1b(input.begin(), input.end(), true);
	const BorelSet bs2(input.begin(), input.end());
	EXPECT_TRUE(bs1a == bs1b);
	EXPECT_FALSE(bs1a == bs2);
	EXPECT_TRUE(bs1a == bs1a);
	EXPECT_TRUE(bs2 == bs2);
	EXPECT_TRUE(bs1a != bs2);
	EXPECT_FALSE(bs1a != bs1b);
}

TEST_F(BorelSetTest,Completion)
{
	std::list<Segment> input1;
	input1.push_back(Segment(0.1, 0.2, false, false));
	input1.push_back(Segment(0.2, 0.25, false, false));
	input1.push_back(Segment(0.4, 0.5, true, false));
	std::list<Segment> input2;
	input2.push_back(Segment(0.1, 0.25, true, true));
	input2.push_back(Segment(0.4, 0.5, true, true));
	const BorelSet bs1(input1.begin(), input1.end());
	const BorelSet bs2(input2.begin(), input2.end());
	EXPECT_EQ(bs2, bs1.complete());
}

TEST_F(BorelSetTest,Union)
{
	std::list<Segment> input1;
	input1.push_back(Segment(0.1, 0.2, false, false));
	input1.push_back(Segment(0.2, 0.25, false, false));
	const BorelSet bs1(input1.begin(), input1.end());
	const BorelSet bs2a(Segment(0.1, 0.2, false, false));
	const BorelSet bs2b(Segment(0.2, 0.22, false, false));
	const BorelSet bs2c(Segment(0.22, 0.25, true, false));
	EXPECT_EQ(bs1, bs2a + bs2b + bs2c);
}

TEST_F(BorelSetTest,Product)
{
	const BorelSet bs1a(Segment(0.1, 0.3, false, false));
	const BorelSet bs1b(Segment(0.15, 0.35, true, true));
	const BorelSet bs2(Segment(0.15, 0.3, true, false));
	EXPECT_EQ(bs2, bs1a * bs1b);
	std::list<Segment> input1;
	input1.push_back(Segment(0.1, 0.2, false, false));
	input1.push_back(Segment(0.2, 0.25, false, false));
	const BorelSet bs3(input1.begin(), input1.end());
	std::list<Segment> input2;
	input2.push_back(Segment(0.15, 0.2, true, false));
	input2.push_back(Segment(0.2, 0.25, false, false));
	const BorelSet bs4(input2.begin(), input2.end());
	EXPECT_EQ(bs4, bs3 * bs2);
}

TEST_F(BorelSetTest,Adjoint)
{
	std::list<Segment> input;
	input.push_back(Segment(-std::numeric_limits<double>::infinity(), 0.1, true, false));
	input.push_back(Segment(0.2, std::numeric_limits<double>::infinity(), true, true));
	const BorelSet expected(input.begin(), input.end());
	const BorelSet bs(Segment(0.1, 0.2, true, false));
	EXPECT_EQ(expected, !bs);
}
