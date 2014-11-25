#include "histogram2d.h"
#include "interpolator2d.h"
#include "random_vector_gaussian.h"
#include "random_vector_levy_symmetric.h"
#include <gtest/gtest.h>
#include <iostream>
#include <Eigen/Core>

using namespace Eigen;

class Histogram2DTest: public testing::Test
{
};

TEST_F(Histogram2DTest,Test)
{
	std::vector<double> x;
	std::vector<double> y;
	// 0 1
	x.push_back(0.1);
	y.push_back(1.2);
	// 0 1
	x.push_back(0.3);
	y.push_back(1.1);
	// 1 1
	x.push_back(1.1);
	y.push_back(1.5);
	// 0 0
	x.push_back(0.4);
	y.push_back(0.41);
	Histogram2D h2d(x, y, 0.0, 2.0, 2, 0.0, 2.0, 2);
	EXPECT_NEAR(1.0, h2d.dx(), 1E-15);
	EXPECT_NEAR(1.0, h2d.dy(), 1E-15);
	ASSERT_EQ(2u, h2d.histogramMatrix().rows());
	ASSERT_EQ(2u, h2d.histogramMatrix().cols());
	EXPECT_EQ(4u, h2d.totalNumber());
	ASSERT_EQ(2u, h2d.xPts().size());
	ASSERT_EQ(2u, h2d.yPts().size());
	for (unsigned int i = 0; i < 2; ++i) {
		EXPECT_NEAR(0.0 + (i+0.5)*h2d.dx(), h2d.xPts()[i], 1E-15);
		EXPECT_NEAR(0.0 + (i+0.5)*h2d.dy(), h2d.yPts()[i], 1E-15);
	}
	EXPECT_NEAR(0.25, h2d.histogramMatrix()(0,0), 1E-15);
	EXPECT_NEAR(0.5, h2d.histogramMatrix()(0,1), 1E-15);
	EXPECT_NEAR(0.0, h2d.histogramMatrix()(1,0), 1E-15);
	EXPECT_NEAR(0.25, h2d.histogramMatrix()(1,1), 1E-15);
	EXPECT_NEAR(1.0, h2d.norm(), 1E-15);
	/*
	Interpolator2D interp(h2d.interpolator());
	EXPECT_NEAR(h2d.xPts().front(), interp.minX(), 1E-15);
	EXPECT_NEAR(h2d.xPts().back(), interp.maxX(), 1E-15);
	EXPECT_NEAR(h2d.yPts().front(), interp.minY(), 1E-15);
	EXPECT_NEAR(h2d.yPts().back(), interp.maxY(), 1E-15);
	for (unsigned int i = 0; i < 2; ++i) {
		for (unsigned int j = 0; j < 2; ++j) {
			EXPECT_NEAR(h2d.histogramMatrix()(i, j), interp(h2d.xPts()[i], h2d.yPts()[j]) * h2d.dx() * h2d.dy(), 1E-15);
		}
	}
	*/
	/*
	for (unsigned int i = 0; i <= 10; ++i) {
		for (unsigned int j = 0; j <= 10; ++j) {
			std::cout << interp(interp.minX() + 0.1*i*(interp.maxX() - interp.minX()), interp.minY() + 0.1*j*(interp.maxY() - interp.minY())) << "\t";
		}
		std::cout << std::endl;
	}
	std::cout << interp.maxX() - interp.minX() << "\t" << interp.maxY() - interp.minY() << "\n";
	*/
	//EXPECT_NEAR(1.0, h2d.interpolator().norm(100, 100, 1), 1E-10);
	//EXPECT_NEAR(1.0, h2d.interpolator(0.5, 2.0, -0.5, 0.5).norm(100, 100, 1), 1E-10);
	/*
	EXPECT_NEAR(0.0, Histogram2D::distance(h2d, 0, 1, 0, 1, h2d, 0, 1, 0, 1, 10, 10, 2, true), 1E-15);
	EXPECT_NEAR(0.0, Histogram2D::distance(h2d, 0, 1, 0, 1, h2d, 0, 1, 0, 1, 10, 10, 2, false), 1E-15);
	EXPECT_NEAR(0.0, Histogram2D::distance(h2d, 0, 2, 0, 3, h2d, 0, 2, 0, 3, 10, 10, 2, true), 1E-15);
	EXPECT_NEAR(0.0, Histogram2D::distance(h2d, 0, 2, 0, 3, h2d, 0, 2, 0, 3, 10, 10, 2, false), 1E-15);
	EXPECT_TRUE(Histogram2D::distance(h2d, 0, 1, 0, 1, h2d, 0, 2, 0, 2, 10, 10, 2, true) > 0);
	EXPECT_TRUE(Histogram2D::distance(h2d, 0, 1, 0, 1, h2d, 0, 2, 0, 2, 10, 10, 2, false) > 0);
	EXPECT_NEAR(2, Histogram2D::distance(h2d, 5, 1, 0, 1, h2d, 0, 1, 0, 1, 100, 100, 1, true), 3E-3);
	EXPECT_NEAR(2, Histogram2D::distance(h2d, 0, 1, 5, 1, h2d, 0, 1, 0, 1, 100, 100, 1, true), 3E-3);
	EXPECT_NEAR(2, Histogram2D::distance(h2d, 0, 1, 0, 1, h2d, 5, 1, 0, 1, 100, 100, 1, true), 3E-3);
	EXPECT_NEAR(2, Histogram2D::distance(h2d, 0, 1, 0, 1, h2d, 0, 1, 5, 1, 100, 100, 1, true), 3E-3);
	*/
}


