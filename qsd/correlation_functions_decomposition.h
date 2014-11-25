#ifndef __QSD_CORRELATION_FUNCTIONS_DECOMPOSITION_H
#define __QSD_CORRELATION_FUNCTIONS_DECOMPOSITION_H

#include <complex>
#include <vector>
#include <math/Jagged2DArray.h>

class CorrelationFunctionsDecomposition
{
public:
	CorrelationFunctionsDecomposition(const rql::math::Jagged2DArray<std::complex<double> >& scales, const rql::math::Jagged2DArray<std::complex<double> >& exponents);
	size_t nbr_sites() const { return m_scales.size(); }
	size_t total_nbr_peaks() const { return m_scales.nbr_elements(); }
	size_t nbr_peaks(size_t site_idx) const { return m_scales.row_size(site_idx); }
	const std::complex<double>& scale(size_t site_idx, size_t peak_idx) const { return m_scales(site_idx, peak_idx); }
	const std::complex<double>& exponent(size_t site_idx, size_t peak_idx) const { return m_exponents(site_idx, peak_idx); }
	std::vector<size_t> peak_counts() const;
private:
	//! alpha_m(tau >= 0) ~= sum_k m_scales(m, k) * exp(tau * m_exponents(m, k))
	rql::math::Jagged2DArray<std::complex<double> > m_scales;
	rql::math::Jagged2DArray<std::complex<double> > m_exponents;
};


#endif // __QSD_CORRELATION_FUNCTIONS_DECOMPOSITION_H
