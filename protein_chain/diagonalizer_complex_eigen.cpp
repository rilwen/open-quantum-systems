#include "diagonalizer_complex_eigen.h"

DiagonalizerComplexEigen::DiagonalizerComplexEigen(unsigned int dim)
	: m_solver(dim)
{
}

void DiagonalizerComplexEigen::diagonalize(const Eigen::MatrixXcd& m, Eigen::MatrixXcd& eigenvectors, Eigen::VectorXcd& eigenvalues)
{
	m_solver.compute(m);
	eigenvectors = m_solver.eigenvectors();
	eigenvalues = m_solver.eigenvalues();
}
