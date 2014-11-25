#include "diagonalizer_realhermitian_eigen.h"

DiagonalizerRealHermitianEigen::DiagonalizerRealHermitianEigen(unsigned int dim)
	: m_dim(dim), m_solver(dim)
{
}

void DiagonalizerRealHermitianEigen::diagonalize(const Eigen::MatrixXd& m, Eigen::MatrixXd& eigenvectors, Eigen::VectorXd& eigenvalues)
{
	m_solver.compute(m);
	eigenvectors = m_solver.eigenvectors();
	eigenvalues = m_solver.eigenvalues();
}
