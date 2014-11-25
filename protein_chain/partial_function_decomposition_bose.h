#ifndef __PARTIAL_FUNCTION_DECOMPOSITION_BOSE_H
#define __PARTIAL_FUNCTION_DECOMPOSITION_BOSE_H

#include <complex>
#include <vector>
#include "core.h"

//! Based on "Hierarchical theory of quantum dissipation: Partial fraction decomposition scheme", Xu J. et al. (2009)
class PartialFunctionDecompositionBose
{
public:
	PROTEIN_CHAIN_API explicit PartialFunctionDecompositionBose(unsigned int order);
	const std::complex<double>& xi(unsigned int n) const
	{
		return m_xi[n];
	}
	PROTEIN_CHAIN_API std::complex<double> evaluate(const std::complex<double>& z) const;
	PROTEIN_CHAIN_API double full_evaluate(double x) const;
	PROTEIN_CHAIN_API_DEBUG std::complex<double> full_evaluate(const std::complex<double>& x) const;
	unsigned int order() const
	{
		return m_order;
	}
private:
	unsigned int m_order;
	std::vector<std::complex<double> > m_xi;
	//std::vector<std::complex<double> > m_sqrt_xi;
};

#endif // __PARTIAL_FUNCTION_DECOMPOSITION_H
