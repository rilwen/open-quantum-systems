#include "correlation_functions_decomposition_factory.h"
#include "correlation_functions_decomposition.h"
#include <protein_chain/correlation_function_lorentzians.h>
#include <protein_chain/correlation_function_decomposable.h>
#include <stdexcept>

namespace CorrelationFunctionsDecompositionFactory
{
	CorrelationFunctionsDecomposition decomposeLorentzians(const std::vector<boost::shared_ptr<const CorrelationFunctionLorentzians> >& alphas)
	{	
		const size_t nbr_sites = alphas.size();
		std::vector<size_t> pc(nbr_sites);
		for (size_t i = 0; i < nbr_sites; ++i) {
			if (!alphas[i])
				throw std::domain_error("Null alpha");
			if (alphas[i]->kBT()==0)
				pc[i] = alphas[i]->gammas().size();
			else
				throw std::domain_error("Only T == 0 supported");
		}
		rql::math::Jagged2DArray<std::complex<double> > scales(pc.begin(), pc.end());
		rql::math::Jagged2DArray<std::complex<double> > exponents(pc.begin(), pc.end());
		for (size_t m = 0; m < nbr_sites; ++m) {
			for (size_t k = 0; k < pc[m]; ++k) {
				exponents(m, k) = std::complex<double>(-alphas[m]->gammas()[k], -alphas[m]->omegas()[k]);
				scales(m, k) = alphas[m]->scales()[k];
			}
		}
		return CorrelationFunctionsDecomposition(scales, exponents);
	}

	CorrelationFunctionsDecomposition decomposeVirtual(const std::vector<boost::shared_ptr<const CorrelationFunctionDecomposable> >& alphas)
	{
		const size_t nbr_sites = alphas.size();
		std::vector<size_t> pc(nbr_sites);
		for (size_t i = 0; i < nbr_sites; ++i) {
			if (!alphas[i])
				throw std::domain_error("Null alpha");
			pc[i] = alphas[i]->nbr_exponents();
		}
		rql::math::Jagged2DArray<std::complex<double> > scales(pc.begin(), pc.end());
		rql::math::Jagged2DArray<std::complex<double> > exponents(pc.begin(), pc.end());
		for (size_t m = 0; m < nbr_sites; ++m) {
			for (size_t k = 0; k < pc[m]; ++k) {
				exponents(m, k) = alphas[m]->exponent(k);
				scales(m, k) = alphas[m]->scale(k).real();
			}
		}
		return CorrelationFunctionsDecomposition(scales, exponents);
	}
}
