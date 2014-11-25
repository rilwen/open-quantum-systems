#ifndef __RQL_MATH_COVARIANCEDECOMPOSITIONCALCULATORITERATIVE_H
#define __RQL_MATH_COVARIANCEDECOMPOSITIONCALCULATORITERATIVE_H

#include "covariance_decomposition_calculator.h"
#include "../MathCore.h"
#include <iostream>
#include <stdexcept>
#include <cassert>
#include <boost/math/special_functions/fpclassify.hpp>
#include <Eigen/Eigenvalues>

namespace rql {
	namespace math {
		namespace stats {
			//! Decomposes covariance matrix iteratively into a transform B such that B_ik = 0 for k > i
			class CovarianceDecompositionCalculatorIterative: public CovarianceDecompositionCalculator
			{
			public:
				RQL_MATH_API CovarianceDecompositionCalculatorIterative();
				RQL_MATH_API void decompose(const Eigen::MatrixXd& covariance, Eigen::MatrixXd& transform) const
				{
					assert( covariance == covariance.adjoint() );
					decompose<Eigen::MatrixXd>(covariance, covariance.rows(), transform);
				}
				template <class M> void decompose(const M& covariance, size_t dimension, Eigen::MatrixXd& transform) const;
				RQL_MATH_API void decompose(const Eigen::MatrixXcd& covariance, Eigen::MatrixXcd& transform) const;
			};

			template <class M> void CovarianceDecompositionCalculatorIterative::decompose(const M& covariance, const size_t dimension, Eigen::MatrixXd& transform) const
			{
				transform.setZero(dimension, dimension);
				if (dimension == 0) {
					return;
				}
				{
					Eigen::ComplexEigenSolver<Eigen::MatrixXcd> diagonalizer;
					Eigen::MatrixXcd covmatrix(dimension, dimension);
					for (size_t r = 0; r < dimension; ++r) {
						for (size_t c = 0; c < dimension; ++c) {
							covmatrix(r, c) = covariance(r, c);
						}
					}
					std::cout << "Covariance hermicity check: " << (covmatrix - covmatrix.adjoint()).norm() << "\n";
					diagonalizer.compute(covmatrix);
					for (size_t i = 0; i < dimension; ++i) {
						std::cout << "Eigenvalue " << i << ": " << diagonalizer.eigenvalues()[i] << "\n";
					}
				}
				assert( covariance(0,0) >= 0 );
				transform(0, 0) = sqrt(covariance(0,0));
				Eigen::MatrixXd covariance_column(dimension, 1);
				// calculate the transform row-by-row, filling in the lower diagonal
				for (size_t i = 1; i < dimension; ++i) {
					covariance_column.setZero();
					// copy covariance values into covariance_column
					for (size_t j = 0; j < i; ++j)
						covariance_column(j, 0) = covariance(j, i);
					// solve the linear equation C_ij = sum_{k=0}^j B_jk B_ik for j < i
					transform.block(i, 0, 1, i) = transform.block(0, 0, i, i).triangularView<Eigen::Lower>().solve(covariance_column.block(0,0,i,1)).transpose();
					// calculate sum_{k=0}^{i-1} (B_ik)^2
					const double tmp1 = transform.block(i, 0, 1, i).cwiseProduct(transform.block(i, 0, 1, i)).sum();
					if (!boost::math::isfinite(tmp1)) {
						throw std::runtime_error("Can't solve for the transform: tmp1 not finite");
					}
					double tmp2 = covariance(i, i) - tmp1;
					if (tmp2 < 0) {
						tmp2 = 0;
					}
					if (!boost::math::isfinite(tmp2)) {
						throw std::runtime_error("Can't solve for the transform: tmp2 not finite");
					}
					if (tmp2 < 0)
						throw std::runtime_error("Can't solve for the transform: tmp2 < 0");
					transform(i, i) = sqrt(tmp2);					
				}
			}
		}
	}
}

#endif // __RQL_MATH_COVARIANCEDECOMPOSITIONCALCULATORITERATIVE_H
