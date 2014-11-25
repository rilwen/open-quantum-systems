#include <gtest/gtest.h>
#include <iostream>
#include "interpolator2d.h"

using namespace Eigen;

class Interpolator2DTest: public testing::Test
{
};

TEST_F(Interpolator2DTest,Flat)
{
	std::vector<double> x(2);
	std::vector<double> y(2);
	MatrixXd z(2,2);
	x[0] = 0;
	x[1] = 1;
	y[0] = 0;
	y[1] = 1;
	const double z0 = 0.5;
	z.fill(z0);
	Interpolator2D intp(x, y, z);
	EXPECT_NEAR(z0, intp(0.1, 0.5), 0.0);
	EXPECT_NEAR(z0, intp(0, 0), 0.0);
	EXPECT_NEAR(z0, intp(0, 1), 0.0);
	EXPECT_NEAR(z0, intp(1, 0), 0.0);
	EXPECT_NEAR(z0, intp(1, 1), 0.0);
}

TEST_F(Interpolator2DTest,SlopeX)
{
	std::vector<double> x(2);
	std::vector<double> y(2);
	MatrixXd z(2,2);
	x[0] = 0;
	x[1] = 1;
	y[0] = 0;
	y[1] = 1;
	const double z0 = 0.5;
	z(0,0) = 0;
	z(0,1) = 0;
	z(1,0) = 1;
	z(1,1) = 1;
	Interpolator2D intp(x, y, z);
	EXPECT_NEAR(0.5, intp(0.5, 0.1), 0.0);
	EXPECT_NEAR(0.5, intp(0.5, 0.7), 0.0);
	EXPECT_NEAR(0.3, intp(0.3, 0.7), 0.0);
	EXPECT_NEAR(0, intp(0, 0), 0.0);
	EXPECT_NEAR(0, intp(0, 1), 0.0);
	EXPECT_NEAR(1, intp(1, 0), 0.0);
	EXPECT_NEAR(1, intp(1, 1), 0.0);
}

TEST_F(Interpolator2DTest,SlopeY)
{
	std::vector<double> x(2);
	std::vector<double> y(2);
	MatrixXd z(2,2);
	x[0] = 0;
	x[1] = 1;
	y[0] = 0;
	y[1] = 1;
	const double z0 = 0.5;
	z(0,0) = 0;
	z(0,1) = 1;
	z(1,0) = 0;
	z(1,1) = 1;
	Interpolator2D intp(x, y, z);
	EXPECT_NEAR(0.5, intp(0.1, 0.5), 0.0);
	EXPECT_NEAR(0.5, intp(0.7, 0.5), 0.0);
	EXPECT_NEAR(0.3, intp(0.7, 0.3), 0.0);
	EXPECT_NEAR(0, intp(0, 0), 0.0);
	EXPECT_NEAR(1, intp(0, 1), 0.0);
	EXPECT_NEAR(0, intp(1, 0), 0.0);
	EXPECT_NEAR(1, intp(1, 1), 0.0);
}

TEST_F(Interpolator2DTest,SlopeXandY)
{
	std::vector<double> x(2);
	std::vector<double> y(2);
	MatrixXd z(2,2);
	x[0] = 0;
	x[1] = 1;
	y[0] = 0;
	y[1] = 1;
	const double z0 = 0.5;
	z(0,0) = 0;
	z(0,1) = 1;
	z(1,0) = 2;
	z(1,1) = 3;
	Interpolator2D intp(x, y, z);
	EXPECT_NEAR(1.5, intp(0.5, 0.5), 0.0);
	EXPECT_NEAR(0.5, intp(0, 0.5), 0.0);
	EXPECT_NEAR(1.4, intp(0.7, 0), 0.0);
	EXPECT_NEAR(1.9, intp(0.7, 0.5), 0.0);
	EXPECT_NEAR(0, intp(0, 0), 0.0);
	EXPECT_NEAR(1, intp(0, 1), 0.0);
	EXPECT_NEAR(2, intp(1, 0), 0.0);
	EXPECT_NEAR(3, intp(1, 1), 0.0);
}

TEST_F(Interpolator2DTest,Multi1)
{
	const unsigned int dim = 4;
	std::vector<double> x(dim);
	std::vector<double> y(dim);
	MatrixXd z(dim,dim);
	const double ax = -0.1;
	const double ay = 0.45;
	const double b = 0.21;
	for (unsigned int i = 0; i < dim; ++i) {
		x[i] = i;
		for (unsigned int j = 0; j < dim; ++j) {
			y[j] = j;
			z(i, j) = b + ax*x[i] + ay*y[j];
		}
	}
	Interpolator2D intp(x, y, z);
	const double cx = 0.4*dim;
	const double cy = 0.1*dim;
	EXPECT_NEAR(b + ax*cx + ay*cy, intp(cx, cy), 1E-15);
	Interpolator2D copy(intp);
	copy *= 0.1;
	EXPECT_NEAR(0.1*intp(cx, cy), copy(cx, cy), 1E-15);
}

TEST_F(Interpolator2DTest,MultiX)
{
	std::vector<double> x(3);
	std::vector<double> y(2);
	MatrixXd z(3,2);
	x[0] = 0;
	x[1] = 1;
	x[2] = 4;
	y[0] = 0;
	y[1] = 1;
	z(0,0) = 0;
	z(1,0) = 0.5;
	z(2,0) = 0.6;
	for (unsigned int i = 0; i < 3; ++i) {
		z(i,1) = z(i,0) + 0.9;
	}
	Interpolator2D intp(x, y, z);
	EXPECT_NEAR(0.55 + 0.1*0.9, intp(2.5, 0.1), 1E-15);
}

TEST_F(Interpolator2DTest,MultiY)
{
	std::vector<double> y(3);
	std::vector<double> x(2);
	MatrixXd z(2,3);
	y[0] = 0;
	y[1] = 1;
	y[2] = 4;
	x[0] = 0;
	x[1] = 1;
	z(0,0) = 0;
	z(0,1) = 0.5;
	z(0,2) = 0.6;
	for (unsigned int i = 0; i < 3; ++i) {
		z(1,i) = z(0,i) + 0.9;
	}
	Interpolator2D intp(x, y, z);
	EXPECT_NEAR(0.55 + 0.1*0.9, intp(0.1,2.5), 1E-15);
}

TEST_F(Interpolator2DTest,ScaleAndShift)
{
	const unsigned int dim = 4;
	std::vector<double> x(dim);
	std::vector<double> y(dim);
	MatrixXd z(dim,dim);
	const double ax = -0.1;
	const double ay = 0.45;
	const double b = 0.21;
	for (unsigned int i = 0; i < dim; ++i) {
		x[i] = i;
		for (unsigned int j = 0; j < dim; ++j) {
			y[j] = j;
			z(i, j) = b + ax*x[i] + ay*y[j];
		}
	}
	Interpolator2D intp1(x, y, z);
	Interpolator2D intp2(intp1);
	EXPECT_NEAR(intp1(0.3*dim,0.41*dim), intp2(0.3*dim,0.41*dim), 1E-20);
	intp2.scale()[0]=2;
	intp2.scale()[1]=0.5;
	EXPECT_NEAR(intp1.minX()/2, intp2.minX(), 1E-15);
	EXPECT_NEAR(intp1.maxX()/2, intp2.maxX(), 1E-15);
	EXPECT_NEAR(intp1.minY()*2, intp2.minY(), 1E-15);
	EXPECT_NEAR(intp1.maxY()*2, intp2.maxY(), 1E-15);
	EXPECT_NEAR(intp1(0.3*dim,0.41*dim), intp2(0.15*dim,0.82*dim), 1E-20);
	intp2.scale().fill(1);
	intp2.shift()[0] = -0.1;
	intp2.shift()[1] = 0.1;
	EXPECT_NEAR(intp1.minX()+0.1, intp2.minX(), 1E-15);
	EXPECT_NEAR(intp1.minY()-0.1, intp2.minY(), 1E-15);
	EXPECT_NEAR(intp1.maxX()+0.1, intp2.maxX(), 1E-15);
	EXPECT_NEAR(intp1.maxY()-0.1, intp2.maxY(), 1E-15);
	EXPECT_NEAR(intp1(0.3*dim,0.41*dim), intp2(0.3*dim + 0.1,0.41 * dim - 0.1), 1E-20);
	intp2.scale()[0]=2;
	intp2.scale()[1]=0.5;
	EXPECT_NEAR(intp1(0.3*dim,0.41*dim), intp2((0.3*dim + 0.1)/2,(0.41 * dim - 0.1)*2), 1E-20);
}


