#include "NormalDistributionTest.h"
#include "math/integr/AdaptiveGK.h"
#include "math/sf/NormalDistribution.h"
#include <iostream>

using namespace rql::math;

// Use sf::erf instead of erf, etc., because C++11 provides its own erf and we do not want a clash.

TEST_F(NormalDistributionTest, CDF)
{
  EXPECT_EQ(0.5, sf::normcdf(0));
	EXPECT_EQ(0, sf::normcdf(-std::numeric_limits<double>::max()));
	EXPECT_EQ(1, sf::normcdf(std::numeric_limits<double>::max()));
	const double x0 = -37.53;
	const double x1 = 8.29;
	EXPECT_NE(0, sf::normcdf(x0));
	EXPECT_NE(1, sf::normcdf(x1));
	for (int i = 0; i <= 100; ++i) {
		const double x = x0 + i*(x1 - x0) / 100.0;
		EXPECT_EQ(sf::normcdf(x), 0.5*sf::erfc(-0.7071067811865475244008443621*x));
	}
	std::cout << sf::normcdf((2.636630902183101-0.06848215242239623)/0.21287763557454142) << std::endl;
}

TEST_F(NormalDistributionTest, Erfc)
{
	const double tolerance1 = 1E-16;
	EXPECT_NEAR(1.999999999999999999999958162, sf::erfc(-7.0), tolerance1);
	EXPECT_NEAR(1.999999999999999978480263288, sf::erfc(-6.0), tolerance1);
	EXPECT_NEAR(1.999999999998462540205571965, sf::erfc(-5.0), tolerance1);
	EXPECT_NEAR(1.999999984582742099719981148, sf::erfc(-4.0), tolerance1);
	EXPECT_NEAR(1.999977909503001414558627224, sf::erfc(-3.0), tolerance1);
	EXPECT_NEAR(1.995322265018952734162069256, sf::erfc(-2.0), tolerance1);
	EXPECT_NEAR(1.842700792949714869341220635, sf::erfc(-1.0), tolerance1);
	EXPECT_NEAR(1.711155633653515131598937835, sf::erfc(-.75), tolerance1);
	EXPECT_NEAR(1.520499877813046537682746654, sf::erfc(-.50), tolerance1);
	EXPECT_NEAR(1.276326390168236932985068268, sf::erfc(-.25), tolerance1);
	EXPECT_EQ(1, sf::erfc(0));
	EXPECT_EQ(2, sf::erfc(-std::numeric_limits<double>::max()));
	EXPECT_EQ(0, sf::erfc(std::numeric_limits<double>::max()));

	EXPECT_NEAR(0.7236736098317630670149317322, sf::erfc(.25), tolerance1);
	EXPECT_NEAR(0.4795001221869534623172533461, sf::erfc(.50), tolerance1);
	EXPECT_NEAR(0.28884436634648486840106216541, sf::erfc(.75), tolerance1);
	EXPECT_NEAR(0.15729920705028513065877936492, sf::erfc(1.0), tolerance1);
	EXPECT_NEAR(0.0046777349810472658379307436330, sf::erfc(2.0), tolerance1);

	// Back to tolerance 1
	EXPECT_NEAR(0.000022090496998585441372776129583, sf::erfc(3.0), tolerance1);
	EXPECT_NEAR(0.000000015417257900280018852159673486, sf::erfc(4.0), tolerance1);
	EXPECT_NEAR(1.5374597944280348501883434853E-12, sf::erfc(5.0), tolerance1);
	EXPECT_NEAR(2.1519736712498913116593350399E-17, sf::erfc(6.0), tolerance1);

	// Tail
	EXPECT_NEAR(4.1838256077794143986140102238E-23, sf::erfc(7.0), 1E-30);
	EXPECT_NEAR(1.12242971729829270799678884432E-29, sf::erfc(8.0), 1E-40);
	EXPECT_NEAR(4.1370317465138102380539034672E-37, sf::erfc(9.0), 1E-50);
	EXPECT_NEAR(2.0884875837625447570007862949E-45, sf::erfc(10.0), 1E-60);

	EXPECT_NE(2, sf::erfc(-5.86));
	EXPECT_NE(0, sf::erfc(26.54));
}

TEST_F(NormalDistributionTest, Erf)
{
	const double tolerance1 = 2E-16;

	EXPECT_NEAR(-0.9999999999999999999999999999887757028270, sf::erf(-8.0), tolerance1);
	EXPECT_NEAR(-0.9999999999999999999999581617439222058560, sf::erf(-7.0), tolerance1);
	EXPECT_NEAR(-0.9999999999999999784802632875010868834066, sf::erf(-6.0), tolerance1);
	EXPECT_NEAR(-0.9999999999984625402055719651498116565146, sf::erf(-5.0), tolerance1);
	EXPECT_NEAR(-0.99999998458274209971998114784032651332, sf::erf(-4.0), tolerance1);
	EXPECT_NEAR(-0.99997790950300141455862722387041706621, sf::erf(-3.0), tolerance1);
	EXPECT_NEAR(-0.99532226501895273416206925636704571137, sf::erf(-2.0), tolerance1);

	EXPECT_NEAR(-0.84270079294971486934122063508, sf::erf(-1.0), 2E-16);
	EXPECT_NEAR(-0.71115563365351513159893783458, sf::erf(-.75), tolerance1);
	EXPECT_NEAR(-0.52049987781304653768274665390, sf::erf(-.50), tolerance1);
	EXPECT_NEAR(-0.27632639016823693298506826776, sf::erf(-.25), 1E-16);
	EXPECT_NEAR(0.84270079294971486934122063508, sf::erf(1.0), 2E-16);
	EXPECT_NEAR(0.71115563365351513159893783458, sf::erf(.75), tolerance1);
	EXPECT_NEAR(0.52049987781304653768274665390, sf::erf(.50), tolerance1);
	EXPECT_NEAR(0.27632639016823693298506826776, sf::erf(.25), 1E-16);

	EXPECT_NEAR(0.9999999999999999999999999999887757028270, sf::erf(8.0), tolerance1);
	EXPECT_NEAR(0.9999999999999999999999581617439222058560, sf::erf(7.0), tolerance1);
	EXPECT_NEAR(0.9999999999999999784802632875010868834066, sf::erf(6.0), tolerance1);
	EXPECT_NEAR(0.9999999999984625402055719651498116565146, sf::erf(5.0), tolerance1);
	EXPECT_NEAR(0.99999998458274209971998114784032651332, sf::erf(4.0), tolerance1);
	EXPECT_NEAR(0.99997790950300141455862722387041706621, sf::erf(3.0), tolerance1);
	EXPECT_NEAR(0.99532226501895273416206925636704571137, sf::erf(2.0), tolerance1);


	EXPECT_EQ(0, sf::erf(0));
	EXPECT_EQ(-1.0, sf::erf(-std::numeric_limits<double>::max()));
	EXPECT_EQ(1.0, sf::erf(std::numeric_limits<double>::max()));

	const double x0 = -5.88;
	const double x1 = 5.88;
	EXPECT_NE(-1, rql::math::sf::erf(x0));
	EXPECT_NE(1, rql::math::sf::erf(x1));
}

TEST_F(NormalDistributionTest, NormsInv)
{
	EXPECT_EQ(-std::numeric_limits<double>::infinity(), sf::normsinv(0));
	EXPECT_EQ(std::numeric_limits<double>::infinity(), sf::normsinv(1));
	EXPECT_EQ(0, sf::normsinv(0.5));
	const double tol = 1E-15;
	std::vector<double> xs;
	xs.push_back(-6);
	xs.push_back(-4.5);
	xs.push_back(-1);
	xs.push_back(0);
	xs.push_back(1);
	xs.push_back(4.5);
	xs.push_back(6);
	for (std::vector<double>::const_iterator i = xs.begin(); i != xs.end(); ++i) {
		const double p = sf::normcdf(*i);
		const double x = sf::normsinv(p);
		const double diff = std::min( std::abs(*i - x), std::abs(p - sf::normcdf(x)));
		EXPECT_NEAR(0, diff, tol);
	}
}

TEST_F(NormalDistributionTest, NormPdf)
{
	rql::math::integr::AdaptiveGK15 igk15(20, 1E-12);
	const double result = igk15.integrate(static_cast<double(*)(double)>(sf::normpdf), -1, 1);
	EXPECT_NEAR(sf::normcdf(1) - sf::normcdf(-1), result, 1E-12);
}
