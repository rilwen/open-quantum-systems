#ifndef __DIAGONALIZER_COMPLEX_H
#define __DIAGONALIZER_COMPLEX_H

#include <Eigen/Core>

class DiagonalizerComplex
{
public:
	virtual void diagonalize(const Eigen::MatrixXcd& m, Eigen::MatrixXcd& eigenvectors, Eigen::VectorXcd& eigenvalues) = 0;
	virtual ~DiagonalizerComplex(){};
};

#endif // __DIAGONALIZER_COMPLEX_H
