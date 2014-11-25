#ifndef __FFTW_COMPLEX_H
#define __FFTW_COMPLEX_H

#include <fftw3.h>
#include <cassert>
#include "utils.h"
#include <Eigen/Core>
#include "core.h"

class FFTWComplex
{
public:
	FFTWComplex(unsigned int size, int sign = FFTW_FORWARD, unsigned int flags = FFTW_MEASURE);
	PROTEIN_CHAIN_API ~FFTWComplex();
	template <class V> void setIn(const V& in);
	PROTEIN_CHAIN_API void setIn(const Eigen::VectorXd& in);
	template <class V> void getOut(V& out) const;
	PROTEIN_CHAIN_API void getOut(unsigned int idx, std::complex<double>& out) const;
	PROTEIN_CHAIN_API void compute();
	template <class V1, class V2> void compute(const V1& in, V2& out);
	unsigned int size() const { return m_size; }
private:
	unsigned int m_size;
	fftw_plan m_plan;
	fftw_complex* m_in;
	fftw_complex* m_out;
};

template <class V> void FFTWComplex::setIn(const V& in)
{
	assert(in.size() == m_size);
	for (unsigned int i = 0; i < m_size; ++i) {
		m_in[i][0] = in[i].real();
		m_in[i][1] = in[i].imag();
	}
}

template <class V> void FFTWComplex::getOut(V& out) const
{
	assert(out.size() == m_size);
	for (unsigned int i = 0; i < m_size; ++i) {
		setReal(out[i], m_out[i][0]);
		setImag(out[i], m_out[i][1]);
	}
}

template <class V1, class V2> inline void FFTWComplex::compute(const V1& in, V2& out)
{
	setIn(in);
	compute();
	getOut(out);
}

#endif // __FFTW_COMPLEX_H
