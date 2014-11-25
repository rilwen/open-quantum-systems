#include "fftw_complex.h"
#include <stdexcept>

FFTWComplex::FFTWComplex(unsigned int size, int sign, unsigned int flags)
	: m_size(size)
{
	m_in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*m_size);
	if (!m_in) {
		// memory allocation failed
		throw std::bad_alloc();
	}
	m_out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*m_size);
	if (!m_out) {
		// memory allocation failed
		throw std::bad_alloc();
	}
	m_plan = fftw_plan_dft_1d(m_size, m_in, m_out, sign, flags);
}

FFTWComplex::~FFTWComplex()
{
	fftw_destroy_plan(m_plan);
	fftw_free(m_in);
	fftw_free(m_out);
}

void FFTWComplex::compute()
{
	fftw_execute(m_plan);
}

void FFTWComplex::setIn(const Eigen::VectorXd& in)
{
	assert(in.size() == m_size);
	for (unsigned int i = 0; i < m_size; ++i) {
		m_in[i][0] = in[i];
		m_in[i][1] = 0.0;
	}
}

void FFTWComplex::getOut(unsigned int idx, std::complex<double>& out) const
{
	assert( idx < m_size );
	setReal(out, m_out[idx][0]);
	setImag(out, m_out[idx][1]);
}
