#include "covariance_decomposition_calculator_one_shot.h"
#include <cassert>
#include <stdexcept>
#include <Eigen/Eigenvalues>

namespace rql {
	namespace math {
		namespace stats {

			void CovarianceDecompositionCalculatorOneShot::decompose(const Eigen::MatrixXd& covariance, Eigen::MatrixXd& transform) const
			{
				assert( covariance == covariance.adjoint() );
				Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver(covariance);
				transform = solver.eigenvectors();
				for (int n = 0; n < covariance.rows(); ++n) {
					const double eval = solver.eigenvalues()[n];
					if (eval >= 0.0)
						transform.col(n) *= sqrt(eval);
					else
						throw std::runtime_error("Negative eigenvalue");
				}
			}

		}
	}
}
