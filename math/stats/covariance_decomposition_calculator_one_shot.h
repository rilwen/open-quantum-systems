#ifndef __RQL_MATH_COVARIANCEDECOMPOSITIONCALCULATORONESHOT_H
#define __RQL_MATH_COVARIANCEDECOMPOSITIONCALCULATORONESHOT_H

#include "covariance_decomposition_calculator.h"
#include "../MathCore.h"

namespace rql {
	namespace math {
		namespace stats {
			class CovarianceDecompositionCalculatorOneShot: public CovarianceDecompositionCalculator
			{
			public:
				RQL_MATH_API void decompose(const Eigen::MatrixXd& covariance, Eigen::MatrixXd& transform) const;
			};
		}
	}
}

#endif // __RQL_MATH_COVARIANCEDECOMPOSITIONCALCULATORONESHOT_H
