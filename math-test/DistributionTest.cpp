#include <vector>
#include <list>
#include "DistributionTest.h"
#include "math/stats/Distribution.h"

using namespace rql::math::stats;

double id(double x)
{
	return x;
}

double id2(double x, const double& f)
{
	return x * f;
}

TEST_F(DistributionTest, ExpectedValue)
{
	std::vector<double> rates(3);
	rates[0] = 0.1;
	rates[1] = 0.4;
	rates[2] = 0.7;
	std::vector<double> rate_probs(3);
	rate_probs[0] = 0.35;
	rate_probs[1] = 0.3;
	rate_probs[2] = 0.35;
	const Distribution<double> dr(rates, rate_probs);
	EXPECT_NEAR(0.4, dr.expected_value(id), 1E-16);
	const double f = 2;
	EXPECT_NEAR(f*0.4, dr.expected_value(id2, f), 1E-15);

	EXPECT_NEAR(0.4, dr.expected_value(), 1E-16);
}

TEST_F(DistributionTest, Construction)
{
	std::vector<int> l1(2);
	l1[0] = 0;
	l1[1] = 1;

	std::vector<double> p1(2);
	p1[0] = 0.8;
	p1[1] = 0.2;

	std::vector<double> p2(2);
	p2[0] = 0.5;
	p2[1] = 0.5;

	const Distribution<int> def_distr;
	EXPECT_EQ(1u, def_distr.size());
	EXPECT_EQ(0, def_distr.value(0));
	EXPECT_EQ(1, def_distr.probability(0));

	const Distribution<int> d1(l1, p1);
	const Distribution<int> d2(p2, 2);
	EXPECT_NEAR(0.8, d1.probability(0), 1E-20);
	EXPECT_NEAR(0.2, d1.probability(1), 1E-20);
	EXPECT_EQ(0, d2.value(0));
	EXPECT_EQ(2, d2.value(1));

	const Distribution<int> single(2);
	EXPECT_EQ(1u, single.size());
	EXPECT_EQ(2, single.value(0));
	EXPECT_EQ(1, single.probability(0));
}

TEST_F(DistributionTest, Convolution)
{
	std::vector<int> l1(2);
	l1[0] = 0;
	l1[1] = 1;

	std::vector<double> p1(2);
	p1[0] = 0.8;
	p1[1] = 0.2;

	std::vector<double> p2(2);
	p2[0] = 0.5;
	p2[1] = 0.5;

	const Distribution<int> d1(l1, p1);
	const Distribution<int> d2(p2, 2);

	std::vector<Distribution<int> > distros;
	distros.push_back(d1);
	distros.push_back(d2);
	Distribution<int> d12 = Distribution<int>::convolute_discrete_positive(distros);
	std::vector<double> v12(4, 0.0);
	Distribution<int>::convolute_discrete_positive(distros, v12);
	EXPECT_EQ(4u, d12.size());
	EXPECT_NEAR(0.8*0.5, d12.probability(0), 1E-16);
	EXPECT_NEAR(0.2*0.5, d12.probability(1), 1E-16);
	EXPECT_NEAR(0.8*0.5, d12.probability(2), 1E-16);
	EXPECT_NEAR(0.2*0.5, d12.probability(3), 1E-16);
	EXPECT_EQ(d12.size(), v12.size());
	for (size_t i = 0; i < d12.size(); ++i) {
		EXPECT_EQ(static_cast<int>(i), d12.value(i));
		EXPECT_EQ(d12.probability(i), v12.at(i));
	}
	Distribution<int> d12copy = d1.convolute_discrete_positive(d2);
	EXPECT_EQ(d12.size(), d12copy.size());
	for (size_t i = 0; i < d12.size(); ++i) {
		EXPECT_EQ(d12.value(i), d12copy.value(i));
		EXPECT_EQ(d12.probability(i), d12copy.probability(i));
	}

	distros.clear();
	distros.push_back(d2);
	distros.push_back(d1);
	d12 = Distribution<int>::convolute_discrete_positive(distros);
	EXPECT_EQ(4u, d12.size());
	EXPECT_NEAR(0.8*0.5, d12.probability(0), 1E-16);
	EXPECT_NEAR(0.2*0.5, d12.probability(1), 1E-16);
	EXPECT_NEAR(0.8*0.5, d12.probability(2), 1E-16);
	EXPECT_NEAR(0.2*0.5, d12.probability(3), 1E-16);

	Distribution<int> single(1);
	distros.clear();
	distros.push_back(d1);
	distros.push_back(single);
	d12 = Distribution<int>::convolute_discrete_positive(distros);
	EXPECT_EQ(d1.size(), d12.size());
	EXPECT_EQ(d1.max() + single.value(0), d12.max());
	for (size_t i = 0; i < d1.size(); ++i) {
		EXPECT_EQ(d1.value(i) + single.value(0), d12.value(i));
		EXPECT_NEAR(d1.probability(i), d12.probability(i), 1E-16);
	}
	distros.clear();
	distros.push_back(single);
	distros.push_back(d1);
	d12 = Distribution<int>::convolute_discrete_positive(distros);
	EXPECT_EQ(d1.size(), d12.size());
	EXPECT_EQ(d1.max() + single.value(0), d12.max());
	for (size_t i = 0; i < d1.size(); ++i) {
		EXPECT_EQ(d1.value(i) + single.value(0), d12.value(i));
		EXPECT_NEAR(d1.probability(i), d12.probability(i), 1E-16);
	}

	std::vector<int> triple_vals(3);
	std::vector<double> triple_probs(3);
	triple_vals[0] = 1;
	triple_vals[1] = 4;
	triple_vals[2] = 5;
	triple_probs[0] = 0.1;
	triple_probs[1] = 0.35;
	triple_probs[2] = 0.55;
	const Distribution<int> d3(triple_vals, triple_probs);
	distros.clear();
	distros.push_back(d3);
	distros.push_back(single);
	v12.resize(1 + d3.max() + single.max());
	std::fill(v12.begin(), v12.end(), 0.0);
	d12 = Distribution<int>::convolute_discrete_positive(distros);
	Distribution<int>::convolute_discrete_positive(distros, v12);
	EXPECT_EQ(1 + d3.max() + single.max(), static_cast<int>(v12.size())); // check the the array hasn't been resized
	EXPECT_EQ(d3.size(), d12.size());
	EXPECT_EQ(d3.max() + single.value(0), d12.max());
	EXPECT_EQ(d3.min() + single.value(0), d12.min());
	EXPECT_NEAR(0.0, v12[0], 1E-16);
	for (size_t i = 0; i < d3.size(); ++i) {
		EXPECT_EQ(d3.value(i) + single.value(0), d12.value(i));
		EXPECT_NEAR(d3.probability(i), d12.probability(i), 1E-16);
		EXPECT_NEAR(d12.probability(i), v12.at(d12.value(i)), 1E-16);
	}
	distros.clear();
	distros.push_back(single);
	distros.push_back(d3);
	std::fill(v12.begin(), v12.end(), 0.0);
	Distribution<int>::convolute_discrete_positive(distros, v12);
	d12 = Distribution<int>::convolute_discrete_positive(distros);
	EXPECT_EQ(d3.size(), d12.size());
	EXPECT_EQ(d3.max() + single.value(0), d12.max());
	EXPECT_EQ(d3.min() + single.value(0), d12.min());
	for (size_t i = 0; i < d3.size(); ++i) {
		EXPECT_EQ(d3.value(i) + single.value(0), d12.value(i));
		EXPECT_NEAR(d3.probability(i), d12.probability(i), 1E-16);
		EXPECT_NEAR(d12.probability(i), v12.at(d12.value(i)), 1E-16);
	}
}

TEST_F(DistributionTest, ConvolutionUnsigned)
{
	std::vector<unsigned int> l1(2);
	l1[0] = 0;
	l1[1] = 1;

	std::vector<double> p1(2);
	p1[0] = 0.8;
	p1[1] = 0.2;

	std::vector<double> p2(2);
	p2[0] = 0.5;
	p2[1] = 0.5;

	const Distribution<unsigned int> d1(l1, p1);
	const Distribution<unsigned int> d2(p2, 2);

	std::vector<Distribution<unsigned int> > distros;
	distros.push_back(d1);
	distros.push_back(d2);
	Distribution<unsigned int> d12 = Distribution<unsigned int>::convolute_discrete_positive(distros);
	std::vector<double> v12(4, 0.0);
	Distribution<unsigned int>::convolute_discrete_positive(distros, v12);
	EXPECT_EQ(4u, d12.size());
	EXPECT_NEAR(0.8*0.5, d12.probability(0), 1E-16);
	EXPECT_NEAR(0.2*0.5, d12.probability(1), 1E-16);
	EXPECT_NEAR(0.8*0.5, d12.probability(2), 1E-16);
	EXPECT_NEAR(0.2*0.5, d12.probability(3), 1E-16);
	EXPECT_EQ(d12.size(), v12.size());
	for (size_t i = 0; i < d12.size(); ++i) {
		EXPECT_EQ(i, d12.value(i));
		EXPECT_EQ(d12.probability(i), v12.at(i));
	}
	Distribution<unsigned int> d12copy = d1.convolute_discrete_positive(d2);
	EXPECT_EQ(d12.size(), d12copy.size());
	for (size_t i = 0; i < d12.size(); ++i) {
		EXPECT_EQ(d12.value(i), d12copy.value(i));
		EXPECT_EQ(d12.probability(i), d12copy.probability(i));
	}

	Distribution<unsigned int> single(1);
	distros.clear();
	distros.push_back(d1);
	distros.push_back(single);
	d12 = Distribution<unsigned int>::convolute_discrete_positive(distros);
	EXPECT_EQ(d1.size(), d12.size());
	EXPECT_EQ(d1.max() + single.value(0), d12.max());
	for (size_t i = 0; i < d1.size(); ++i) {
		EXPECT_EQ(d1.value(i) + single.value(0), d12.value(i));
		EXPECT_NEAR(d1.probability(i), d12.probability(i), 1E-16);
	}
	distros.clear();
	distros.push_back(single);
	distros.push_back(d1);
	d12 = Distribution<unsigned int>::convolute_discrete_positive(distros);
	EXPECT_EQ(d1.size(), d12.size());
	EXPECT_EQ(d1.max() + single.value(0), d12.max());
	for (size_t i = 0; i < d1.size(); ++i) {
		EXPECT_EQ(d1.value(i) + single.value(0), d12.value(i));
		EXPECT_NEAR(d1.probability(i), d12.probability(i), 1E-16);
	}

	std::vector<unsigned int> triple_vals(3);
	std::vector<double> triple_probs(3);
	triple_vals[0] = 1;
	triple_vals[1] = 4;
	triple_vals[2] = 5;
	triple_probs[0] = 0.1;
	triple_probs[1] = 0.35;
	triple_probs[2] = 0.55;
	const Distribution<unsigned int> d3(triple_vals, triple_probs);
	distros.clear();
	distros.push_back(d3);
	distros.push_back(single);
	v12.resize(1 + d3.max() + single.max());
	std::fill(v12.begin(), v12.end(), 0.0);
	d12 = Distribution<unsigned int>::convolute_discrete_positive(distros);
	Distribution<unsigned int>::convolute_discrete_positive(distros, v12);
	EXPECT_EQ(1 + d3.max() + single.max(), v12.size()); // check the the array hasn't been resized
	EXPECT_EQ(d3.size(), d12.size());
	EXPECT_EQ(d3.max() + single.value(0), d12.max());
	EXPECT_EQ(d3.min() + single.value(0), d12.min());
	for (size_t i = 0; i < d12.min(); ++i) {
	    EXPECT_NEAR(0.0, v12[i], 1E-16);
	}
	for (size_t i = 0; i < d3.size(); ++i) {
		EXPECT_EQ(d3.value(i) + single.value(0), d12.value(i));
		EXPECT_NEAR(d3.probability(i), d12.probability(i), 1E-16);
		EXPECT_NEAR(d12.probability(i), v12.at(d12.value(i)), 1E-16);
	}
	distros.clear();
	distros.push_back(single);
	distros.push_back(d3);
	std::fill(v12.begin(), v12.end(), 0.0);
	d12 = Distribution<unsigned int>::convolute_discrete_positive(distros);
	Distribution<unsigned int>::convolute_discrete_positive(distros, v12);
	EXPECT_EQ(d3.size(), d12.size());
	EXPECT_EQ(d3.max() + single.value(0), d12.max());
	EXPECT_EQ(d3.min() + single.value(0), d12.min());
	for (size_t i = 0; i < d3.size(); ++i) {
		EXPECT_EQ(d3.value(i) + single.value(0), d12.value(i));
		EXPECT_NEAR(d3.probability(i), d12.probability(i), 1E-16);
		EXPECT_NEAR(d12.probability(i), v12.at(d12.value(i)), 1E-16);
	}
}

TEST_F(DistributionTest, Copy)
{
	std::vector<unsigned int> l1(2);
	l1[0] = 0;
	l1[1] = 1;

	std::vector<double> p1(2);
	p1[0] = 0.8;
	p1[1] = 0.2;

	const Distribution<unsigned int> d1(l1, p1);

	std::list<unsigned int> copy_l1;
	std::list<double> copy_p1;
	d1.copy_values(copy_l1);
	d1.copy_probabilities(copy_p1);
	EXPECT_EQ(2u, copy_l1.size());
	EXPECT_EQ(2u, copy_p1.size());
	std::list<unsigned int>::iterator il1 = copy_l1.begin();
	EXPECT_EQ(l1[0], *il1);
	++il1;
	EXPECT_EQ(l1[1], *il1);
	std::list<double>::iterator ip1 = copy_p1.begin();
	EXPECT_EQ(p1[0], *ip1);
	++ip1;
	EXPECT_EQ(p1[1], *ip1);

	std::vector<unsigned int> vc_l1(2);
	std::vector<double> vc_p1(2);
	std::copy(d1.values_begin(), d1.values_end(), vc_l1.begin());
	std::copy(d1.probabilities_begin(), d1.probabilities_end(), vc_p1.begin());
	for (size_t i = 0; i < d1.size(); ++i) {
	    EXPECT_EQ(l1[i], vc_l1[i]);
	    EXPECT_EQ(p1[i], vc_p1[i]);
	}
}
