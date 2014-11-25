#include "BrentTest.h"
#include "math/solvers/Brent.h"
#include "math/sf/NormalDistribution.h"
#include "math/MathUtils.h"

using namespace rql::math::solvers;
using namespace rql::math::sf;

int counter;

double xinvcoshx(double x)
{
	++counter;
	return x != 0 ? x/cosh(x) : 0;
}

TEST_F(BrentTest, XInvCoshX)
{
	const double tol = 1E-18;
	counter = 0;
	EXPECT_NEAR(0, find_root(xinvcoshx, -1, 5, tol, tol), tol);
	std::cout << "Used " << counter << " function evaluations" << std::endl;
//	std::cout << "Proceeded by bisection " << maths::root_finder_stats.bisection_cnt << " times" << std::endl;
//	std::cout << "Proceeded by inverse quadratic interpolation " << maths::root_finder_stats.inverse_quadratic_cnt << " times" << std::endl;
//	std::cout << "Proceeded by linear interpolation " << maths::root_finder_stats.linear_cnt << " times" << std::endl;
}

double my_atan(double x)
{
	++counter;
	return atan(x);
}

TEST_F(BrentTest, ATan)
{
	const double tol = 1E-18;
	counter = 0;
	EXPECT_NEAR(0, find_root(my_atan, -1, 5, tol, tol), tol);
	std::cout << "Used " << counter << " function evaluations" << std::endl;
//	std::cout << "Proceeded by bisection " << maths::root_finder_stats.bisection_cnt << " times" << std::endl;
//	std::cout << "Proceeded by inverse quadratic interpolation " << maths::root_finder_stats.inverse_quadratic_cnt << " times" << std::endl;
//	std::cout << "Proceeded by linear interpolation " << maths::root_finder_stats.linear_cnt << " times" << std::endl;
}

double polynomial(double x)
{
	++counter;
	return (x - 1.0001)*(x - 1)*(x-4)*(x-10);
}

TEST_F(BrentTest, Polynomial)
{
	const double tol = 1E-18;
	counter = 0;
	const double x0 = find_root(polynomial, -7, 7, tol, tol);
	std::cout << "Used " << counter << " function evaluations" << std::endl;
	EXPECT_NEAR(0, polynomial(x0), effective_tolerance(1, tol));
//	std::cout << "Proceeded by bisection " << maths::root_finder_stats.bisection_cnt << " times" << std::endl;
//	std::cout << "Proceeded by inverse quadratic interpolation " << maths::root_finder_stats.inverse_quadratic_cnt << " times" << std::endl;
//	std::cout << "Proceeded by linear interpolation " << maths::root_finder_stats.linear_cnt << " times" << std::endl;
//	counter = 0;
//	const double brent_x0 = brent(polynomial, -7, 7, tol, tol);
//	std::cout << "Used " << counter << " function evaluations" << std::endl;
//	EXPECT_NEAR(0, polynomial(brent_x0), effective_tolerance(1, tol));
}

double lines(double x)
{
	++counter;
	if (x < -1) {
		return -1;
	} else if (x < 1) {
		return x;
	}
	else return 1;
}

TEST_F(BrentTest, Lines)
{
	const double tol = 1E-18;
	counter = 0;
	EXPECT_NEAR(0, find_root(lines, -10, 5, tol, tol), tol);
	std::cout << "Used " << counter << " function evaluations" << std::endl;
//	std::cout << "Proceeded by bisection " << maths::root_finder_stats.bisection_cnt << " times" << std::endl;
//	std::cout << "Proceeded by inverse quadratic interpolation " << maths::root_finder_stats.inverse_quadratic_cnt << " times" << std::endl;
//	std::cout << "Proceeded by linear interpolation " << maths::root_finder_stats.linear_cnt << " times" << std::endl;
}

struct NormCdfInverter
{
	NormCdfInverter(double new_p)
	{
		p = new_p;
		counter = 0;
	}
	double p;
	mutable int counter;
	double operator()(double x) const
	{
		++counter;
		return normcdf(x) - p;
	}
};

TEST_F(BrentTest, NormCDF)
{
	const double tol = 1E-18;
	const double x0 = -3.5;
	NormCdfInverter inverter(normcdf(x0));
	EXPECT_NEAR(x0, find_root(inverter, -11, 10, tol, tol), effective_tolerance(x0, tol));
	std::cout << "Used " << inverter.counter << " function evaluations" << std::endl;
//	std::cout << "Proceeded by bisection " << maths::root_finder_stats.bisection_cnt << " times" << std::endl;
//	std::cout << "Proceeded by inverse quadratic interpolation " << maths::root_finder_stats.inverse_quadratic_cnt << " times" << std::endl;
//	std::cout << "Proceeded by linear interpolation " << maths::root_finder_stats.linear_cnt << " times" << std::endl;
//	std::cout << "Rejected interpolation " << maths::root_finder_stats.rejection_cnt << " times" << std::endl;
}

// See http://orion.math.iastate.edu/burkardt/f_src/testzero/testzero.html

double p02(double x)
{
	++counter;
//	std::cout << "Evaluating at x == " << std::setprecision(16) << x << std::endl;
	return 2*x - exp(-x);
}

TEST_F(BrentTest, P02)
{
	const double tol = 1E-18;
	counter = 0;
	double xa_tolx;
	double xb_tolx;
	double xa_toly;
	double xb_toly;
//	std::cout << "TOLX" << std::endl;
	const double root_tolx = find_root(p02, -10, 100, tol, 0, xa_tolx, xb_tolx);
	std::cout << "Used " << counter << " function evaluations" << std::endl;
//	std::cout << "Proceeded by bisection " << maths::root_finder_stats.bisection_cnt << " times" << std::endl;
//	std::cout << "Proceeded by inverse quadratic interpolation " << maths::root_finder_stats.inverse_quadratic_cnt << " times" << std::endl;
//	std::cout << "Proceeded by linear interpolation " << maths::root_finder_stats.linear_cnt << " times" << std::endl;
	EXPECT_NEAR(0, xb_tolx - xa_tolx, effective_tolerance(root_tolx, tol));
//	std::cout << "xb - xa == " << xb_tolx - xa_tolx << std::endl;
//	std::cout << "f(root) == " << p02(root_tolx) << std::endl;
//	std::cout << "f(xa) == " << p02(xa_tolx) << std::endl;
//	std::cout << "f(xb) == " << p02(xb_tolx) << std::endl;
//	counter = 0;
//	const double brent_root_tolx = brent(p02, -10, 100, tol, 0);
//	std::cout << "Used " << counter << " function evaluations for Brent" << std::endl;

	counter = 0;
//	std::cout << "TOLY" << std::endl;
	const double root_toly = find_root(p02, -10, 100, 0, tol, xa_toly, xb_toly);
	std::cout << "Used " << counter << " function evaluations" << std::endl;
//	std::cout << "Proceeded by bisection " << maths::root_finder_stats.bisection_cnt << " times" << std::endl;
//	std::cout << "Proceeded by inverse quadratic interpolation " << maths::root_finder_stats.inverse_quadratic_cnt << " times" << std::endl;
//	std::cout << "Proceeded by linear interpolation " << maths::root_finder_stats.linear_cnt << " times" << std::endl;
	EXPECT_NEAR(0, p02(root_toly), tol);
//	std::cout << "xb - xa == " << xb_toly - xa_toly << std::endl;
//	std::cout << "f(root) == " << p02(root_toly) << std::endl;
//	std::cout << "f(xa) == " << p02(xa_toly) << std::endl;
//	std::cout << "f(xb) == " << p02(xb_toly) << std::endl;

}

double pinhead(double x)
{
	static const double epsilon = 0.00001E+00;
	++counter;
	return ( 16.0E+00 - x*x*x*x ) / ( 16.0E+00 * x*x*x*x + epsilon );
}

TEST_F(BrentTest, Pinhead)
{
	const double tol = 1E-18;
	counter = 0;
	double xa_tolx, xb_tolx;
	double xa_toly, xb_toly;
	const double root_tolx = find_root(pinhead, 0, 10, tol, 0, xa_tolx, xb_tolx);
	std::cout << "Used " << counter << " function evaluations" << std::endl;
	EXPECT_NEAR(0, xb_tolx - xa_tolx, effective_tolerance(root_tolx, tol));
	EXPECT_NEAR(2, root_tolx, effective_tolerance(root_tolx, tol));
	counter = 0;
	const double root_toly = find_root(pinhead, 0, 10, 0, tol, xa_toly, xb_toly);
	std::cout << "Used " << counter << " function evaluations" << std::endl;
	EXPECT_NEAR(0, pinhead(root_toly), tol);

//	counter = 0;
//	const double root_brent_tolx = brent(pinhead, 0, 10, tol, 0);
//	std::cout << "Used " << counter << " function evaluations for Brent" << std::endl;
//	std::cout << "Function value at Brent root: " << pinhead(root_brent_tolx) << std::endl;
//	counter = 0;
//	const double root_brent_toly = brent(pinhead, 0, 10, 0, tol);
//	std::cout << "Used " << counter << " function evaluations for Brent" << std::endl;
//	std::cout << "Function value at Brent root: " << pinhead(root_brent_toly) << std::endl;
}

double repeller(double x)
{
	++counter;
	return 20.0E+00 * x / ( 100.0E+00 * x * x + 1.0E+00 );
}

TEST_F(BrentTest, Repeller)
{
	const double tol = 1E-18;
	counter = 0;
	double xa_tolx;
	double xb_tolx;
	double xa_toly;
	double xb_toly;
	const double root_tolx = find_root(repeller, -11, 8, tol, 0, xa_tolx, xb_tolx);
	std::cout << "Used " << counter << " function evaluations" << std::endl;
	EXPECT_NEAR(0, xb_tolx - xa_tolx, effective_tolerance(root_tolx, tol));
	EXPECT_NEAR(0, root_tolx, effective_tolerance(root_tolx, tol));
//	counter = 0;
//	const double root_brent_tolx = brent(repeller, -11, 8, tol, 0);
//	std::cout << "Used " << counter << " function evaluations for Brent" << std::endl;
	counter = 0;
	const double root_toly = find_root(repeller, -11, 8, 0, tol, xa_toly, xb_toly);
	std::cout << "Used " << counter << " function evaluations" << std::endl;
	EXPECT_NEAR(0, repeller(root_toly), tol);
//	counter = 0;
//	const double root_brent_toly = brent(repeller, -11, 8, 0, tol);
//	std::cout << "Used " << counter << " function evaluations for Brent" << std::endl;
}

double camel(double x)
{
	++counter;
	return 1.0E+00 / ( pow( x - 0.3E+00, 2 ) + 0.01E+00 )
	+ 1.0E+00 / ( pow( x - 0.9E+00, 2 ) + 0.04E+00 ) + 2.0E+00 * x - 5.2E+00;
}

TEST_F(BrentTest, Camel)
{
	const double tol = 1E-18;
	counter = 0;
	double xa_tolx, xb_tolx;
//	double xa_toly;
//	double xb_toly;
	const double root_tolx = find_root(camel, -11, 10, tol, 0, xa_tolx, xb_tolx);
	std::cout << "Used " << counter << " function evaluations" << std::endl;
//	std::cout << "Proceeded by bisection " << maths::root_finder_stats.bisection_cnt << " times" << std::endl;
//	std::cout << "Proceeded by inverse quadratic interpolation " << maths::root_finder_stats.inverse_quadratic_cnt << " times" << std::endl;
//	std::cout << "Proceeded by linear interpolation " << maths::root_finder_stats.linear_cnt << " times" << std::endl;
//	std::cout << "Rejected interpolation " << maths::root_finder_stats.rejection_cnt << " times" << std::endl;
//	std::cout << "Function value at our root: " << camel(root_tolx) << std::endl;
	EXPECT_NEAR(0, xb_tolx - xa_tolx, effective_tolerance(root_tolx, tol));
	counter = 0;
//	const double root_toly = find_root(camel, -11, 10, 0, tol, xa_toly, xb_toly);
	std::cout << "Used " << counter << " function evaluations" << std::endl;
//	std::cout << "Proceeded by bisection " << maths::root_finder_stats.bisection_cnt << " times" << std::endl;
//	std::cout << "Proceeded by inverse quadratic interpolation " << maths::root_finder_stats.inverse_quadratic_cnt << " times" << std::endl;
//	std::cout << "Proceeded by linear interpolation " << maths::root_finder_stats.linear_cnt << " times" << std::endl;
//	std::cout << "Rejected interpolation " << maths::root_finder_stats.rejection_cnt << " times" << std::endl;
//	std::cout << "Function value at our root: " << camel(root_toly) << std::endl;
	EXPECT_NEAR(0, xb_tolx - xa_tolx, effective_tolerance(root_tolx, 0));
//	counter = 0;
//	const double root_brent_tolx = brent(camel, -10, 10, tol, 0);
//	std::cout << "Used " << counter << " function evaluations for Brent" << std::endl;
//	std::cout << "Function value at Brent root: " << camel(root_brent_tolx) << std::endl;
//	counter = 0;
//	const double root_brent_toly = brent(camel, -10, 10, 0, tol);
//	std::cout << "Used " << counter << " function evaluations for Brent" << std::endl;
//	std::cout << "Function value at Brent root: " << camel(root_brent_toly) << std::endl;
//	EXPECT_NEAR(0, camel(root_toly), tol); // effective tolerance kicks in first
}

double lazy_boy(double x)
{
	++counter;
	return 0.00000000001E+00 * (x - 100);
}

TEST_F(BrentTest, LazyBoy)
{
	const double tol = 1E-18;
	counter = 0;
	double xa_tolx, xb_tolx;
	double xa_toly, xb_toly;
	const double root_tolx = find_root(lazy_boy, -12000000000000.0E+00, 10000000000000.0E+00, tol, 0, xa_tolx, xb_tolx);
	std::cout << "Used " << counter << " function evaluations" << std::endl;
	EXPECT_NEAR(0, xb_tolx - xa_tolx, effective_tolerance(root_tolx, tol));
	EXPECT_NEAR(100, root_tolx, effective_tolerance(root_tolx, tol));
	counter = 0;
	const double root_toly = find_root(lazy_boy, -12000000000000.0E+00, 10000000000000.0E+00, 0, tol, xa_toly, xb_toly);
	std::cout << "Used " << counter << " function evaluations" << std::endl;
	EXPECT_NEAR(0, lazy_boy(root_toly), tol);

//	counter = 0;
//	const double root_brent_tolx = brent(lazy_boy, -10000000000000.0E+00, 10000000000000.0E+00, tol, 0);
//	std::cout << "Used " << counter << " function evaluations for Brent" << std::endl;
//	std::cout << "Function value at Brent root: " << lazy_boy(root_brent_tolx) << std::endl;
//	counter = 0;
//	const double root_brent_toly = brent(lazy_boy, -10000000000000.0E+00, 10000000000000.0E+00, 0, tol);
//	std::cout << "Used " << counter << " function evaluations for Brent" << std::endl;
//	std::cout << "Function value at Brent root: " << lazy_boy(root_brent_toly) << std::endl;
}

double vertigo(double x)
{
	++counter;
	const double y = x - 1E6;
	const int sign = y > 0 ? 1 : -1;
	return 1E12*sign*sqrt(std::abs(y));
}

TEST_F(BrentTest, Vertigo)
{
	const double tol = 1E-18;
	counter = 0;
	double xa_tolx, xb_tolx;
	double xa_toly, xb_toly;
	const double root_tolx = find_root(vertigo, -2E6, 3E6, tol, 0, xa_tolx, xb_tolx);
	std::cout << "Used " << counter << " function evaluations" << std::endl;
	EXPECT_NEAR(0, xb_tolx - xa_tolx, effective_tolerance(root_tolx, tol));
	EXPECT_NEAR(1E6, root_tolx, effective_tolerance(root_tolx, tol));
	counter = 0;
	find_root(vertigo, -2E6, 3E6, 0, tol, xa_toly, xb_toly);
	std::cout << "Used " << counter << " function evaluations" << std::endl;
	EXPECT_NEAR(0, xb_toly - xa_toly, effective_tolerance(root_tolx, 0));
}

//double extrap1(double xa, double fa, double xb, double fb)
//{
//	return (fa*xb - fb*xa) / (fa - fb);
//}
//
//double extrap2(double xa, double fa, double xb, double fb)
//{
//	double inv_slope = (xb - xa) / (fb - fa);
//	return xa - inv_slope * fa;
//}
//
//double extrap3(double xa, double fa, double xb, double fb)
//{
//	double slope = (fb - fa) / (xb - xa);
//	return xa - fa / slope;
//}
//
//TEST_F(BrentTest, LinearExtrapolation)
//{
//	const double x = 0.5E9;
//	const double f = 1E-14;
//	const double xc1 = extrap1(-x, f, x, -f);
//	const double xc2 = extrap2(-x, f, x, -f);
//	const double xc3 = extrap3(-x, f, x, -f);
//	std::cout << "Error xc1 == " << xc1 << std::endl;
//	std::cout << "Error xc2 == " << xc2 << std::endl;
//	std::cout << "Error xc3 == " << xc3 << std::endl;
//}
