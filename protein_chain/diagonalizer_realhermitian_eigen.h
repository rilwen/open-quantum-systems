#ifndef __DIAGONALIZER_REALHERMITIAN_EIGEN_H
#define __DIAGONALIZER_REALHERMITIAN_EIGEN_H

#include <Eigen/Core>
#include <Eigen/Eigenvalues>
#include "diagonalizer_realhermitian.h"
#include "core.h"

class DiagonalizerRealHermitianEigen: public DiagonalizerRealHermitian
{
public:
	PROTEIN_CHAIN_API DiagonalizerRealHermitianEigen(unsigned int dim);
	PROTEIN_CHAIN_API void diagonalize(const Eigen::MatrixXd& m, Eigen::MatrixXd& eigenvectors, Eigen::VectorXd& eigenvalues);
	PROTEIN_CHAIN_API_DEBUG unsigned int dim() const { return m_dim; }
private:
	unsigned int m_dim;
	Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd>	m_solver;
 };

#endif
