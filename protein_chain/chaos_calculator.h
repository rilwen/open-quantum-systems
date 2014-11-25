#ifndef __CHAOS_CALCULATOR_H
#define __CHAOS_CALCULATOR_H

#include "fftw_complex.h"
#include <Eigen/Core>
#include <math/average.h>
#include <vector>
#include "core.h"

class ChaosCalculator
{
public:
	//! @param size Number of eigenvalues
	PROTEIN_CHAIN_API ChaosCalculator(unsigned int size);
	template <class V> void calculate(const V& eigenvalues);
	const Eigen::VectorXd& spectrum() const { return m_spectrum; }
	const std::vector<rql::math::Average<double> >& averages() const { return m_averages; }
private:
	FFTWComplex m_fft;
	Eigen::VectorXd m_deltas;
	Eigen::VectorXd m_spectrum;
	std::vector<rql::math::Average<double> > m_averages;
};

template <class V> void ChaosCalculator::calculate(const V& eigenvalues) 
{
	if (eigenvalues.size() != m_fft.size() ) throw std::domain_error("eigenvalues.size() != m_fft.size()+1");
	for (unsigned int i = 1; i < m_fft.size(); ++i) {
		m_deltas[i] = eigenvalues[i] - eigenvalues[i-1];
		if (m_deltas[i] < 0) throw std::runtime_error("eigenvalues are not sorted");
	}
	const double s_avg = m_deltas.sum()/(m_deltas.size()-1);
	for (unsigned int i = 1; i < m_fft.size(); ++i) {
		m_deltas[i] -= s_avg;
	}
	for (unsigned int n = 2; n < m_fft.size(); ++n) {
		m_deltas[n] += m_deltas[n-1];
	}

	m_fft.setIn(m_deltas);
	m_fft.compute();
	
	std::complex<double> z; 
	for (unsigned int k = 0; k < m_fft.size(); ++k) {
		m_fft.getOut(k, z);
		m_spectrum[k] = std::abs(z)*std::abs(z)/m_fft.size();
		m_averages[k].update(m_spectrum[k]);
	}

}

#endif // __CHAOS_CALCULATOR_H
