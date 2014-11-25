#ifndef __DIAGONALIZER_REALHERMITIAN_H
#define __DIAGONALIZER_REALHERMITIAN_H

#include <Eigen/Core>

class DiagonalizerRealHermitian
{
public:
	virtual unsigned int dim() const = 0;
	virtual void diagonalize(const Eigen::MatrixXd& m, Eigen::MatrixXd& eigenvectors, Eigen::VectorXd& eigenvalues) = 0;
	virtual ~DiagonalizerRealHermitian(){};
};

#endif
