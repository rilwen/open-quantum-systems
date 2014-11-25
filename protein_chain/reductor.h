#ifndef __REDUCTOR_HPP
#define __REDUCTOR_HPP

#include <Eigen/Core>
#include "core.h"

//! Calculates reduced density matrices
class Reductor
{
public:
	PROTEIN_CHAIN_API Reductor();
	PROTEIN_CHAIN_API Reductor(size_t dimA, size_t dimB);
	size_t dimA() const { return m_dim_A; }
	size_t dimB() const { return m_dim_B; }
	
	//! Reduce |state><state|
	PROTEIN_CHAIN_API void reduce(const Eigen::VectorXcd& state, Eigen::MatrixXcd& rho) const;
	
	//! Reduce |left><right|
	PROTEIN_CHAIN_API void reduce(const Eigen::VectorXcd& left, const Eigen::VectorXcd& right, Eigen::MatrixXcd& rho) const;
private:
	size_t m_dim_A;
	size_t m_dim_B;
	size_t m_basis_dim;
};

#endif //  __REDUCTOR_HPP
