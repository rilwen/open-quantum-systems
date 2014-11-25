#include "hamiltonian.h"
#include <stdexcept>
#include <Eigen/Eigenvalues>
#include <iostream>
#include "evolver_taylor_expansion.h"
#include <cmath>
#include "utils.h"

using namespace Eigen;

Hamiltonian::Hamiltonian(bool isHermitian, const Eigen::MatrixXcd& matrix)
	: m_dim(matrix.rows()), m_isHermitian(isHermitian), m_matrix(matrix), m_eigenvalues(m_dim), m_eigenvectors(m_dim,m_dim)
{
	if (m_dim != matrix.cols()) throw std::domain_error("Matrix must be square");
	m_canDiverge = false;
	if (m_isHermitian) {
		m_isNormal = true;
		SelfAdjointEigenSolver<MatrixXcd> diagonalizer(m_matrix);
		m_eigenvalues = diagonalizer.eigenvalues().cast<std::complex<double> >();
		m_eigenvectors = diagonalizer.eigenvectors();
	} else {
		ComplexEigenSolver<MatrixXcd> diagonalizer(m_matrix);
		m_eigenvalues = diagonalizer.eigenvalues();
		m_eigenvectors = diagonalizer.eigenvectors();
		for (unsigned int i = 0; i < m_dim; ++i) {
			if (m_eigenvalues[i].imag() > 1e-14)  m_canDiverge = true;			
		}
		MatrixXcd::ConjugateReturnType mconj = m_matrix.conjugate();
		const double comnorm = (m_matrix*mconj - mconj*m_matrix).norm();
		m_isNormal = comnorm < 1E-12;
//		std::cout << "comnorm = " << comnorm << std::endl;
	}
	double max_abs_eval = 0;
	for (unsigned int i = 0; i < m_dim; ++i) {
		max_abs_eval = std::max(max_abs_eval, std::abs(m_eigenvalues[i].real()));
		max_abs_eval = std::max(max_abs_eval, std::abs(m_eigenvalues[i].imag()));
	}
	const double time_step = max_abs_eval > 0 ? 1.0/(10*max_abs_eval) : 1;
	m_evolver = boost::shared_ptr<Evolver>(new EvolverTaylorExpansion<4u>(m_matrix, time_step));
}

boost::shared_ptr<Hamiltonian> Hamiltonian::buildChain(const std::vector<double>& epsilons, double J, double kappa)
{
	MatrixXcd m(epsilons.size(), epsilons.size());
	m *= 0.0;

	for (unsigned int i = 0; i < epsilons.size(); ++i) {
		m(i,i) = epsilons[i];
		if (i>0) {
			m(i-1,i) = m(i,i-1) = J;
		}
	}
	setImag(m(epsilons.size() - 1, epsilons.size() - 1), kappa);
	const bool isHermitian = kappa == 0;

	return boost::shared_ptr<Hamiltonian> (new Hamiltonian(isHermitian, m));
}

boost::shared_ptr<Hamiltonian> Hamiltonian::buildDipoles(const std::vector<double>& epsilons, double J, double kappa)
{
	//! kappa is ignored for now
	const bool isHermitian = true;

	MatrixXcd m(epsilons.size(), epsilons.size());
	m *= 0.0;

	for (unsigned int i = 0; i < epsilons.size(); ++i) {
		m(i,i) = epsilons[i];
		for (unsigned int j = 0; j < epsilons.size(); ++j) {
			if (j != i) {
				m(i, j) = J / pow(std::abs(static_cast<double>(i) - j), 3);
			}
		}
	}

	return boost::shared_ptr<Hamiltonian> (new Hamiltonian(isHermitian, m));
}

boost::shared_ptr<Hamiltonian> Hamiltonian::buildRing(const std::vector<double>& epsilons, double J, double kappa)
{
	MatrixXcd m(epsilons.size(), epsilons.size());
	m *= 0.0;

	for (unsigned int i = 0; i < epsilons.size(); ++i) {
		m(i,i) = epsilons[i];
		if (i>0) {
			m(i-1,i) = m(i,i-1) = J;
			m(i-1, epsilons.size()-1) = m(epsilons.size()-1, i-1) = J;
		}
	}
	m(epsilons.size()-2,0) = m(0,epsilons.size()-2) = J;

	setImag(m(epsilons.size() - 1, epsilons.size() - 1), kappa);
	const bool isHermitian = kappa == 0;

	return boost::shared_ptr<Hamiltonian> (new Hamiltonian(isHermitian, m));
}
