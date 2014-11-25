#include <gtest/gtest.h>
#include "eigenstate_properties.h"
#include <Eigen/Core>
#include <cmath>

using namespace Eigen;

class EigenstatePropertiesTest: public testing::Test
{
};

static const double pi = 3.141592653589793;

TEST_F(EigenstatePropertiesTest, PR)
{
	VectorXd v(100);
	for (unsigned int i = 0; i < 100; ++i) {
		v[i] = 0.1;
	}
	EXPECT_NEAR(100.0, EigenstateProperties::pr(v), 1E-10);
	for (unsigned int i = 0; i < 100; ++i) {
		v[i] = sin( pi*34*(i + 1)/101.0 ) * sqrt(2/101.0);
	}
	EXPECT_NEAR(101.0*2/3.0, EigenstateProperties::pr(v), 1E-10);
}

TEST_F(EigenstatePropertiesTest, Mu)
{
	VectorXd v(100);
	for (unsigned int i = 0; i < 100; ++i) {
		v[i] = 0.1;
	}
	EXPECT_NEAR(10.0, EigenstateProperties::mu(v), 1E-10);		
	for (unsigned int i = 1; i < 100; ++i) {
		v[i] = -v[i-1];
	}
	EXPECT_NEAR(0.0, EigenstateProperties::mu(v), 1E-10);		
}

TEST_F(EigenstatePropertiesTest, L1Norm)
{
	VectorXd v(100);
	for (unsigned int i = 0; i < 100; ++i) {
		v[i] = 0.1;
	}
	EXPECT_NEAR(10.0, EigenstateProperties::l1norm(v), 1E-10);		
	for (unsigned int i = 1; i < 100; ++i) {
		v[i] = -v[i-1];
	}
	EXPECT_NEAR(10.0, EigenstateProperties::l1norm(v), 1E-10);		
}
