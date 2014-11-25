#ifndef __RQL_MATH_COVARIANCEDECOMPOSITIONCALCULATOR_H
#define __RQL_MATH_COVARIANCEDECOMPOSITIONCALCULATOR_H

#include <Eigen/Core>

namespace rql {
	namespace math {
		namespace stats {
			class CovarianceDecompositionCalculator
			{
			public:
				virtual ~CovarianceDecompositionCalculator() {}
				//! For given N x N symmetric covariance matrix C, calculate such B that C_ij = sum_{k=1}^N B_ik B_jk
				//! Consequently, for independent N(0,1) variables x_k, y_i = sum_k B_ik x_k have cov(y_i, y_j) = C_ij
				virtual void decompose(const Eigen::MatrixXd& covariance, Eigen::MatrixXd& transform) const = 0;
			};
		}
	}
}

#endif // __RQL_MATH_COVARIANCEDECOMPOSITIONCALCULATOR_H
