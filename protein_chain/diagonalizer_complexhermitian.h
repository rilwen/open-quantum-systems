#ifndef __DIAGONALIZER_COMPLEXHERMITIAN_H
#define __DIAGONALIZER_COMPLEXHERMITIAN_H

#include <Eigen/Core>

class DiagonalizerComplexHermitian
{
public:
	virtual void diagonalize(const Eigen::MatrixXcd& m, Eigen::MatrixXcd& eigenvectors, Eigen::VectorXd& eigenvalues) = 0;
	virtual ~DiagonalizerComplexHermitian(){};
};

#endif // __DIAGONALIZER_COMPLEXHERMITIAN_H
