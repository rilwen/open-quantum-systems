#include "cfa_utils.h"
#include <protein_chain/spectral_density.h>
#include <protein_chain/correlation_function_lorentzians.h>
#include <cmath>
#include <boost/make_shared.hpp>

namespace cfa {
	const std::complex<double> IMAGINARY(0,1);

	void discretise_correlation_function(const std::vector<boost::shared_ptr<const SpectralDensity> >& J, const double w0, const double w1
		, const size_t nbr_bath_modes_per_site, Eigen::VectorXd& omega, Eigen::MatrixXcd& g)
	{
		const size_t nbr_sites = J.size();
		const size_t nbr_bath_modes = nbr_sites * nbr_bath_modes_per_site;
		const double dw = (w1 - w0) / nbr_bath_modes_per_site;
		omega.setZero(nbr_bath_modes);
		g.setZero(nbr_bath_modes, nbr_sites);
		for (size_t m = 0; m < nbr_sites; ++m) {
			const size_t k0 = m*nbr_bath_modes_per_site;
			for (size_t dk = 0; dk < nbr_bath_modes_per_site; ++dk) {
				const double w = w0 + (0.5 + dk) * dw;
				const size_t k = k0 + dk;
				const double spectral = J[m]->spectral_density(w);
				g(k, m) = sqrt(spectral * dw);
				omega[k] = w;
			}
		}
	}

	void discretise_correlation_function(boost::shared_ptr<const SpectralDensity> J, double w0, double w1, size_t nbr_sites, size_t nbr_bath_modes_per_site
		, Eigen::VectorXd& omega, Eigen::MatrixXcd& g)
	{
		std::vector<boost::shared_ptr<const SpectralDensity> > Jvec(nbr_sites, J);
		discretise_correlation_function(Jvec, w0, w1, nbr_bath_modes_per_site, omega, g);
	}

	boost::shared_ptr<CorrelationFunctionLorentzians> discretise_correlation_function(boost::shared_ptr<const SpectralDensity> J, double w0, double w1, size_t nbr_bath_modes_per_site)
	{
		const double dw = (w1 - w0) / nbr_bath_modes_per_site;
		Eigen::VectorXd omega(nbr_bath_modes_per_site);
		Eigen::VectorXd hwhm(nbr_bath_modes_per_site);
		Eigen::VectorXd height(nbr_bath_modes_per_site);
		for (size_t i = 0; i < nbr_bath_modes_per_site; ++i) {
			const double w = w0 + (0.5 + i) * dw;
			omega[i] = w;
			const double curr_hwhm = dw;
			hwhm[i] = curr_hwhm;
			const double scale = J->spectral_density(w) * dw;
			height[i] = CorrelationFunctionLorentzians::height(scale, curr_hwhm);
		}
		return boost::make_shared<CorrelationFunctionLorentzians>(omega, hwhm, height);
	}
}
