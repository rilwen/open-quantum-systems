#include "random_vector_gaussian.h"
#include <math/stats/covariance.h>
#include <boost/math/distributions/normal.hpp>
#include <cassert>
#include <cmath>
#include <limits>
#include <Eigen/SVD>

using namespace Eigen;
using namespace rql::math::stats;

RandomVectorGaussian::RandomVectorGaussian(unsigned int dim, double sigma, size_t batch_size, bool antithetic_sampling)
	: m_dim(dim), m_engine(17), m_uniform_01(), m_uniform_rng(m_engine, m_uniform_01), m_sigma(sigma), m_batch(batch_size * dim), m_batch_iter(m_batch.end()), m_negate(false), m_antithetic_sampling(antithetic_sampling)
{
}

static const double p_low = 0.02425;
static const double p_high = 1 - p_low;
static const boost::math::normal_distribution<double> standard_normal(0, 1);

// Implementation of the algorithm described in http://home.online.no/~pjacklam/notes/invnorm/
double RandomVectorGaussian::normsinv(double p)
{
	assert( p >= 0 );
	assert( p <= 1 );
	int sign = 1;
	if (p > 0.5) {
		p = 1 - p;
		sign = -1;
	}
	if (p == 0) {
		return -sign * std::numeric_limits<double>::infinity();
	}

	static const double a[] = {-3.969683028665376e+01, 2.209460984245205e+02,
		-2.759285104469687e+02, 1.383577518672690e+02,
		-3.066479806614716e+01, 2.506628277459239e+00};
	static const double b[] = {-5.447609879822406e+01, 1.615858368580409e+02,
		-1.556989798598866e+02, 6.680131188771972e+01,
		-1.328068155288572e+01};
	static const double c[] = {-7.784894002430293e-03, -3.223964580411365e-01,
		-2.400758277161838e+00, -2.549732539343734e+00,
		4.374664141464968e+00, 2.938163982698783e+00};
	static const double d[] = {7.784695709041462e-03, 3.224671290700398e-01,
		2.445134137142996e+00, 3.754408661907416e+00};
	double x;
	if (p < p_low) {
		const double q = sqrt(-2*log(p));
		x = (((((c[0]*q+c[1])*q+c[2])*q+c[3])*q+c[4])*q+c[5]) /
			((((d[0]*q+d[1])*q+d[2])*q+d[3])*q+1);
	} else {
		const double q = p - 0.5;
		const double r = q*q;
		x = (((((a[0]*r+a[1])*r+a[2])*r+a[3])*r+a[4])*r+a[5])*q /
			(((((b[0]*r+b[1])*r+b[2])*r+b[3])*r+b[4])*r+1);
	}
	const double e = boost::math::cdf(standard_normal, x) - p;
	if (x > -37) {
		const double u = e * 2.506628274631000502415765285 * exp(0.5*x*x);
		return sign * (x - u/(1 + 0.5*x*u));
	} else {
		return x;
	}
}

template <class V> void RandomVectorGaussian::drawImpl(V& vec)
{
	assert( vec.size() == m_dim );
	if (m_batch_iter == m_batch.end()) {
		draw_batch();
		m_batch_iter = m_batch.begin();
	}
	for (size_t i = 0; i < m_dim; ++i) {
		assert( m_batch_iter != m_batch.end() );
		vec[i] = *m_batch_iter*m_sigma;
		++m_batch_iter;
	}
}

void RandomVectorGaussian::draw(std::vector<double>& vec)
{
	drawImpl(vec);
}

void RandomVectorGaussian::draw(Eigen::VectorXd& vec)
{
	drawImpl(vec);
}

void RandomVectorGaussian::draw_batch()
{
	if (m_negate) {
		for (std::vector<double>::iterator it = m_batch.begin(); it != m_batch.end(); ++it) {			
			*it *= -1;
		}
		m_negate = false;
	} else {
		for (std::vector<double>::iterator it = m_batch.begin(); it != m_batch.end(); ++it) {
			*it = normsinv(m_uniform_rng());
		}
		if (m_antithetic_sampling)
			m_negate = true;
		//decorrelate();
	}
}

void RandomVectorGaussian::decorrelate()
{
	Covariance covcalc(m_dim);
	std::vector<double>::iterator i1, i2;
	i1 = m_batch.begin();
	i2 = i1;
	while (i1 != m_batch.end()) {
		i2 += m_dim;
		covcalc.update(i1, i2);
		i1 = i2;
	}
	MatrixXd cov(covcalc.covariance());
	JacobiSVD<MatrixXd> svd(cov, ComputeFullU | ComputeFullV);
	MatrixXd s(MatrixXd::Zero(m_dim,m_dim));
	for (unsigned int i = 0; i < m_dim; ++i) {
		s(i, i) = 1.0 / sqrt(svd.singularValues()[i]);
	}
	s = s*svd.matrixV().conjugate();
	const MatrixXd transform(svd.matrixU()*s);

	std::vector<double> decorred(m_dim);
	i1 = m_batch.begin();
	while (i1 != m_batch.end()) {
		for (unsigned int k = 0; k < m_dim; ++k) {
			i2 = i1;
			double sum = 0;
			for (unsigned int l = 0; l < m_dim; ++l) {
				//assert( i2 != m_batch.end() );
				sum += transform(l,k) * (*i2);
				++i2;
			}
			decorred[k] = sum;
		}
		i2 = i1;
		i1 += m_dim;
		std::copy(decorred.begin(), decorred.end(), i2);
	}
}

void RandomVectorGaussian::set_seed(unsigned int seed)
{
	m_engine.seed(seed);
}

const size_t RandomVectorGaussian::DEFAULT_BATCH_SIZE = 5000;

