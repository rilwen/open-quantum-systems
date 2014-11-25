#ifndef __CFA_UTILS_H
#define __CFA_UTILS_H

#include <cassert>
#include <complex>
#include <iostream>
#include <vector>
#include <boost/shared_ptr.hpp>
#include <Eigen/Core>
#include "core.h"

class CorrelationFunctionLorentzians;
class SpectralDensity;

namespace cfa {
	inline size_t exciton_index(size_t i, size_t j) 
	{ 
		assert(i >= j);
		return (i*(i + 1))/2 + j; 
	}

	template <class M> bool is_hermitian(const M& matrix, double rel_eps = 1E-6)
	{
		const double m = matrix.norm();
		const double d = (matrix - matrix.adjoint()).norm();
		if (d > rel_eps*m) {
			std::cerr << "NONHERMITIAN: " << matrix << std::endl;
			return false; // you can set a breakpoint here
		} else {
			return true;
		}
	};

	//! Convert A into (A + A^dagger)/2
	template <class M> void make_hermitian(M& matrix)
	{
		const int rsize = matrix.rows();
		//const int csize = matrix.cols();
		for (int r = 0; r < rsize; ++r) {
			matrix(r, r) = matrix(r, r).real();
			for (int c = 0; c < r; ++c) {
				std::complex<double> tmp = 0.5*(matrix(r, c) + conj(matrix(c, r)));
				matrix(r, c) = tmp;
				matrix(c, r) = conj(tmp);
			}
		}
	}

	//! A --> A - A^dagger
	template <class M> void subtract_adjoint_in_place(M& matrix, const size_t dim) 
	{
		for (size_t i = 0; i < dim; ++i) {
			matrix(i, i) = std::complex<double>(0.0, 2*matrix(i, i).imag());
			for (size_t j = 0; j < i; ++j) {
				matrix(i, j) -= conj(matrix(j, i));
				matrix(j, i) = -conj(matrix(i, j));
			}
		}
	}

	inline bool is_zero(const std::complex<double>& z)
	{
		return z.real() == 0 && z.imag() == 0;
	}

	CFA_API void discretise_correlation_function(const std::vector<boost::shared_ptr<const SpectralDensity> >& J, double w0, double w1, size_t nbr_bath_modes_per_site, Eigen::VectorXd& omega, Eigen::MatrixXcd& g);
	CFA_API void discretise_correlation_function(boost::shared_ptr<const SpectralDensity> J, double w0, double w1, size_t nbr_sites, size_t nbr_bath_modes_per_site, Eigen::VectorXd& omega, Eigen::MatrixXcd& g);
	CFA_API boost::shared_ptr<CorrelationFunctionLorentzians> discretise_correlation_function(boost::shared_ptr<const SpectralDensity> J, double w0, double w1, size_t nbr_bath_modes_per_site);

	CFA_API extern const std::complex<double> IMAGINARY;
}

#endif // __CFA_UTILS_H
