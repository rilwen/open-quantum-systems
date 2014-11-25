#include "reductor.h"
#include <cassert>
#include <stdexcept>
#include <iostream>

Reductor::Reductor()
	: m_dim_A(1), m_dim_B(1), m_basis_dim(1)
{
}

Reductor::Reductor(size_t dimA, size_t dimB)
	: m_dim_A(dimA), m_dim_B(dimB), m_basis_dim(dimA*dimB)
{
	assert(dimA);
	assert(dimB);
}

void Reductor::reduce(const Eigen::VectorXcd& state, Eigen::MatrixXcd& rho) const
{
	reduce(state, state, rho);
}

void Reductor::reduce(const Eigen::VectorXcd& left, const Eigen::VectorXcd& right, Eigen::MatrixXcd& rho) const {
	if (left.size() != static_cast<int>(m_basis_dim)) {
		std::cerr << "state size = " << left.size() << ", m_basis_dim = " << m_basis_dim << std::endl;
		throw std::domain_error("Wrong left state dimension");
	}
	if (right.size() != static_cast<int>(m_basis_dim)) {
		std::cerr << "state size = " << right.size() << ", m_basis_dim = " << m_basis_dim << std::endl;
		throw std::domain_error("Wrong right state dimension");
	}
	rho.setZero(m_dim_A, m_dim_A);
	for (size_t i1 = 0; i1 < m_dim_A; ++i1) {
		for (size_t i2 = 0; i2 < m_dim_A; ++i2) {
			for (size_t j = 0; j < m_dim_B; ++j) {
				rho(i1, i2) += conj(left[i1*m_dim_B + j])*right[i2*m_dim_B + j];
			}
		}
	}
}

