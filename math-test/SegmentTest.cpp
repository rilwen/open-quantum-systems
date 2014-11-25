#include "SegmentTest.h"
#include "math/Segment.h"

using namespace rql::math;

TEST_F(SegmentTest,Comparison)
{
	Segment s1(0.1, 0.2, true, false);
	Segment s2(0.1, 0.2, true, false);
	EXPECT_EQ(s1, s2);
	EXPECT_TRUE(s1 == s2);
	EXPECT_FALSE(s1 != s2);
	EXPECT_LE(s1, s2);
	EXPECT_FALSE(s1 < s2);
	Segment s3(0.1, 0.21, true, false);
	EXPECT_NE(s1, s3);
	EXPECT_TRUE(s1 != s3);
	EXPECT_FALSE(s1 == s3);
	EXPECT_LT(s1, s3);
	EXPECT_LE(s1, s3);
	Segment s4(0.05, 0.2, true, true);
	EXPECT_NE(s1, s4);
	EXPECT_LT(s4, s1);
	Segment s5(0.1, 0.2, true, true);
	EXPECT_NE(s1, s5);
	EXPECT_LE(s1, s5);
	EXPECT_LT(s1, s5);
	Segment s6(0.1, 0.2, false, false);
	EXPECT_NE(s1, s6);
	EXPECT_LE(s1, s6);
	EXPECT_LT(s1, s6);
}

TEST_F(SegmentTest,IsAbove)
{
	const Segment s1(0.1, 0.2, true, false);
	const Segment s2a(0.2, 0.25, true, false);
	const Segment s2b(0.2, 0.25, false, false);
	const Segment s2c(0.15, 0.25, false, false);
	EXPECT_TRUE(s1.is_below(s2a));
	EXPECT_TRUE(s1.is_below(s2b));
	EXPECT_FALSE(s1.is_below(s2c));
	EXPECT_FALSE(s2a.is_below(s1));
	EXPECT_FALSE(s2c.is_below(s1));
}

TEST_F(SegmentTest,MinDistance)
{
	const Segment s1a(0.1, 0.2, true, false);
	EXPECT_EQ(0.0, s1a.min_distance(s1a));
	const Segment s1b(0.1, 0.2, false, true);
	EXPECT_EQ(0.0, s1a.min_distance(s1b));
	const Segment s2a(0.2, 0.3, false, true);
	EXPECT_EQ(0.0, s1a.min_distance(s2a));
	const Segment s2b(0.4, 0.6, false, true);
	EXPECT_EQ(0.2, s1a.min_distance(s2b));
	EXPECT_EQ(0.2, s2b.min_distance(s1a));
	const Segment s3(0.15, 0.25, true, true);
	EXPECT_EQ(0.0, s1a.min_distance(s3));
}

TEST_F(SegmentTest,CanBeMergedWith)
{
	const Segment s1a(0.1, 0.2, true, false);
	const Segment s1b(0.16, 0.21, false, false);
	EXPECT_TRUE(s1a.can_be_merged_with(s1a));
	EXPECT_TRUE(s1a.can_be_merged_with(s1b));
	EXPECT_TRUE(s1b.can_be_merged_with(s1b));
	EXPECT_TRUE(s1b.can_be_merged_with(s1a));
	const Segment s2a(0.2, 0.3, true, true);
	const Segment s2b(0.2, 0.3, false, true);
	EXPECT_TRUE(s1a.can_be_merged_with(s2a));
	EXPECT_TRUE(s2a.can_be_merged_with(s1a));
	EXPECT_FALSE(s1a.can_be_merged_with(s2b));
	EXPECT_FALSE(s2b.can_be_merged_with(s1a));
	EXPECT_TRUE(s1a.can_be_merged_with(s2b, true));
	EXPECT_TRUE(s2b.can_be_merged_with(s1a, true));
	EXPECT_TRUE(s1b.can_be_merged_with(s2b));
	EXPECT_TRUE(s2b.can_be_merged_with(s1b));
	const Segment s3(0.25, 0.5, true, true);
	EXPECT_FALSE(s1a.can_be_merged_with(s3));

	EXPECT_TRUE(Segment::EMPTY.can_be_merged_with(Segment::EMPTY));
	EXPECT_TRUE(s1a.can_be_merged_with(Segment::EMPTY));
	EXPECT_TRUE(Segment::EMPTY.can_be_merged_with(s1a));
}

TEST_F(SegmentTest,MergeWith)
{
	const Segment s1a(0.1, 0.2, true, false);
	EXPECT_EQ(s1a, s1a.merge_with(Segment::EMPTY));
	EXPECT_EQ(s1a, Segment::EMPTY.merge_with(s1a));
	EXPECT_EQ(Segment::EMPTY, Segment::EMPTY.merge_with(Segment::EMPTY));
	EXPECT_EQ(s1a, s1a.merge_with(s1a));
	const Segment s1b(0.2, 0.3, true, true);
	EXPECT_EQ(Segment(0.1, 0.3, true, true), s1a.merge_with(s1b));
	EXPECT_EQ(Segment(0.1, 0.3, true, true), s1b.merge_with(s1a));
	const Segment s1c(0.2, 0.3, false, true);
	EXPECT_EQ(Segment(0.1, 0.3, true, true), s1a.merge_with(s1c, true));
	EXPECT_EQ(Segment(0.1, 0.3, true, true), s1c.merge_with(s1a, true));
	const Segment s2a(0.15, 0.2, false, true);
	EXPECT_EQ(Segment(0.1,0.2,true,true), s1a.merge_with(s2a));
	EXPECT_EQ(Segment(0.1,0.2,true,true), s2a.merge_with(s1a));
	const Segment s2b(0.15, 0.25, false, false);
	EXPECT_EQ(Segment(0.1,0.25,true,false), s1a.merge_with(s2b));
	EXPECT_EQ(Segment(0.1,0.25,true,false), s2b.merge_with(s1a));
}

TEST_F(SegmentTest,Complete)
{
	EXPECT_EQ(Segment::EMPTY, Segment::EMPTY.complete());
	const Segment s1a(0.1, 0.2, true, true);
	EXPECT_EQ(s1a, s1a.complete());
	const Segment s1b(0.1, 0.2, false, false);
	EXPECT_EQ(s1a, s1b.complete());
}