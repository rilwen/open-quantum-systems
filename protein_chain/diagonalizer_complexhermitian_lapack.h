#ifndef __DIAGONALIZER_COMPLEXHERMITIAN_LAPACK_H
#define __DIAGONALIZER_COMPLEXHERMITIAN_LAPACK_H

#include "diagonalizer_complexhermitian.h"

class DiagonalizerComplexHermitianLapack: public DiagonalizerComplexHermitian
{
public:
	void diagonalize(const Eigen::MatrixXcd& m, Eigen::MatrixXcd& eigenvectors, Eigen::VectorXd& eigenvalues);
};

#endif // __DIAGONALIZER_COMPLEXHERMITIAN_LAPACK_H
