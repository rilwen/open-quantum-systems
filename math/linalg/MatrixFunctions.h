#ifndef __RQL_MATH_MATRIX_FUNCTIONS_H
#define __RQL_MATH_MATRIX_FUNCTIONS_H

#include <Eigen/Core>
#include <Eigen/Eigenvalues>
#include <cassert>

namespace rql {
	namespace math {
		namespace linalg {
			namespace MatrixFunctions
			{
				namespace {
					template <class M, class F> void hermitian_matrix_function_impl(const M& A, const F function, M& result);
				}

				template <class F> void hermitian_matrix_function(const Eigen::MatrixXcd& A, const F function, Eigen::MatrixXcd& result)
				{
					hermitian_matrix_function_impl(A, function, result);
				}

				template <class F> void hermitian_matrix_function(const Eigen::MatrixXd& A, const F function, Eigen::MatrixXd& result)
				{
					hermitian_matrix_function_impl(A, function, result);
				}

				namespace {
					template <class M, class F> void hermitian_matrix_function_impl(const M& A, const F function, M& result)
					{
						assert(A.rows() == A.cols());
						const int dimension = A.rows();						
						Eigen::SelfAdjointEigenSolver<M> diagonaliser(A);
						result.setZero(dimension, dimension);
						for (int i = 0; i < dimension; ++i) {
							const double result_eigenvalue = function(diagonaliser.eigenvalues()[i]);
							result += result_eigenvalue * diagonaliser.eigenvectors().col(i) * diagonaliser.eigenvectors().col(i).adjoint();
						}
					}
				}

				inline double id(double x)
				{
					return x;
				}
			}
		}
	}
}

#endif // __RQL_MATHS_MATRIX_FUNCTIONS_H
