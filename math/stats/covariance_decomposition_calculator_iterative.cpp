#include "covariance_decomposition_calculator_iterative.h"
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <Eigen/Dense>

namespace rql {
	namespace math {
		namespace stats {
			CovarianceDecompositionCalculatorIterative::CovarianceDecompositionCalculatorIterative()
			{
			}

			void CovarianceDecompositionCalculatorIterative::decompose(const Eigen::MatrixXcd& covariance, Eigen::MatrixXcd& transform) const
			{
				assert( covariance == covariance.adjoint() );
				transform.setZero(covariance.rows(), covariance.rows());
				if (covariance.rows() == 0)
					return;
				// Solve for covariance_kl = sum_j=0^min(k,l) \cc{transform_kj} transform_lj
				for (int k = 0; k < covariance.rows(); ++k) {
					for (int l = 0; l < k; ++l) {
						assert(l < k);
						std::complex<double> tmp(0.0);
						for (int j = 0; j < l; ++j) {
							tmp += conj(transform(k, j)) * transform(l, j);
						}
						transform(k, l) = (covariance(k, l) - tmp) / transform(l, l);
					}
					std::complex<double> tmp(0.0);
					for (int j = 0; j < k; ++j) {
						if (!boost::math::isfinite(transform(k, j))) {
							throw std::runtime_error("Can't solve for complex transform");
						}
						tmp += conj(transform(k, j)) * transform(k, j);
					}
					transform(k, k) = sqrt(covariance(k, k) - tmp);
				}
			}

			//void CovarianceDecompositionCalculatorIterative::decompose(const Eigen::MatrixXd& covariance, Eigen::MatrixXd& transform) const
			//{
			//	assert( covariance == covariance.adjoint() );
			//	transform.setZero(covariance.rows(), covariance.cols());
			//	if (covariance.rows() == 0)
			//		return;
			//	assert( covariance(0,0) >= 0 );
			//	transform(0, 0) = sqrt(covariance(0,0));
			//	// calculate the transform row-by-row, filling in the lower diagonal
			//	for (int i = 1; i < covariance.rows(); ++i) {
			//		// solve the linear equation C_ij = sum_{k=0}^j B_jk B_ik for j < i
			//		transform.block(i, 0, 1, i) = transform.topLeftCorner(i,i).triangularView<Eigen::Lower>().solve(covariance.block(0,i,i,1)).transpose();
			//		// calculate sum_{k=0}^{i-1} (B_ik)^2
			//		const double tmp1 = transform.block(i, 0, 1, i).cwiseProduct(transform.block(i, 0, 1, i)).sum();
			//		if (!boost::math::isfinite(tmp1)) {
			//			throw std::runtime_error("Can't solve for the transform");
			//		}
			//		const double tmp2 = std::max(covariance(i, i) - tmp1, 0.0); // avoiding negative cov matrix values
			//		if (!boost::math::isfinite(tmp1)) {
			//			throw std::runtime_error("Can't solve for the transform");
			//		}
			//		if (tmp2 < 0)
			//			throw std::runtime_error("Can't solve for the transform");
			//		transform(i, i) = sqrt(tmp2);
			//	}
			//}			
		}
	}
}