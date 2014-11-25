#ifndef __DIAGONALIZER_COMPLEX_EIGEN_H
#define __DIAGONALIZER_COMPLEX_EIGEN_H

#include "diagonalizer_complex.h"
#include <Eigen/Eigenvalues>
#include "core.h"

class DiagonalizerComplexEigen: public DiagonalizerComplex
{
public:
	PROTEIN_CHAIN_API DiagonalizerComplexEigen(unsigned int dim);
	PROTEIN_CHAIN_API void diagonalize(const Eigen::MatrixXcd& m, Eigen::MatrixXcd& eigenvectors, Eigen::VectorXcd& eigenvalues);
private:
	Eigen::ComplexEigenSolver<Eigen::MatrixXcd>	m_solver;
};

#endif // __DIAGONALIZER_COMPLEX_EIGEN_H
