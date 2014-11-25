#include "correlation_functions_decomposition.h"
#include <stdexcept>

CorrelationFunctionsDecomposition::CorrelationFunctionsDecomposition(const rql::math::Jagged2DArray<std::complex<double> >& scales, const rql::math::Jagged2DArray<std::complex<double> >& exponents)
	: m_scales(scales), m_exponents(exponents)
{
	const size_t nr = scales.size();
	if (exponents.size() != nr) {
		throw std::domain_error("Row nbr mismatch");
	}
	for (size_t i = 0; i < nr; ++i) {
		if (scales.row_size(i) != exponents.row_size(i)) {
			throw std::domain_error("Col nbr mismatch");
		}
	}
}

std::vector<size_t> CorrelationFunctionsDecomposition::peak_counts() const
{
	std::vector<size_t> pc(nbr_sites());
	for (size_t m = 0; m < pc.size(); ++m)
		pc[m] = nbr_peaks(m);
	return pc;
}
