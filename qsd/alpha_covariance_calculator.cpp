#include "alpha_covariance_calculator.h"
#include <protein_chain/CorrelationFunction.h>
#include <stdexcept>

AlphaCovarianceCalculator::AlphaCovarianceCalculator(const std::vector<boost::shared_ptr<const CorrelationFunction> >& alpha, const size_t nbrTimePts, double d)
	: m_d(d), m_alpha_cov(alpha.size(), nbrTimePts)
{
	for (int m = 0; m < m_alpha_cov.rows(); ++m) {
		for (size_t k = 0; k < nbrTimePts; ++k) {
			m_alpha_cov(m, k) = (*alpha[m])(k*d);
		}
		m_alpha_cov(m, 0) = m_alpha_cov(m, 0).real();
	}
}

AlphaCovarianceCalculator::AlphaCovarianceCalculator(boost::shared_ptr<const CorrelationFunction> alpha, const size_t dim, const size_t nbrTimePts, double d)
	: m_d(d), m_alpha_cov(dim, nbrTimePts)
{
	for (size_t k = 0; k < nbrTimePts; ++k) {
		m_alpha_cov.col(k).fill((*alpha)(k*d));
	}
	for (int m = 0; m < m_alpha_cov.rows(); ++m) {
		m_alpha_cov(m, 0) = m_alpha_cov(m, 0).real();
	}
}

std::complex<double> AlphaCovarianceCalculator::covariance(size_t m, size_t k, size_t l) const
{
	assert( m < static_cast<size_t>(m_alpha_cov.rows()));
	assert( k < static_cast<size_t>(m_alpha_cov.cols()));
	assert( l < static_cast<size_t>(m_alpha_cov.cols()));
	if (k >= l)
		return m_alpha_cov(m, k - l);
	else
		return conj(m_alpha_cov(m, l - k));
}

double AlphaCovarianceCalculator::covariance_RR_II(size_t m, size_t k, size_t l) const
{
	const std::complex<double> cov(covariance(m, k, l));
	return 0.5*cov.real();
}

double AlphaCovarianceCalculator::covariance_RI(size_t m, size_t k, size_t l) const
{
	const std::complex<double> cov(covariance(m, k, l));
	return 0.5*cov.imag();
}

SingleAlphaCovarianceCalculator::SingleAlphaCovarianceCalculator(const AlphaCovarianceCalculator& acc, size_t m)
	: m_acc(acc), m_m(m)
{
}

const SingleAlphaCovarianceCalculator AlphaCovarianceCalculator::single_alpha_calculator(size_t m) const
{
	return SingleAlphaCovarianceCalculator(*this, m);
}
