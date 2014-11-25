#ifndef __HAMILTONIAN_H
#define __HAMILTONIAN_H

#include <Eigen/Core>
#include <boost/shared_ptr.hpp>
#include <vector>
#include "core.h"

class Evolver;

class Hamiltonian
{
public:
	Hamiltonian(bool isHermitian, const Eigen::MatrixXcd& matrix);
	unsigned int dim() const {return m_dim;}
	bool isHermitian() const {return m_isHermitian;}
	bool isNormal() const {return m_isNormal;}
	const Eigen::MatrixXcd& matrix() const {return m_matrix;}
	const Eigen::VectorXcd& eigenvalues() const {return m_eigenvalues;}
	const Eigen::MatrixXcd& eigenvectors() const {return m_eigenvectors;}
	boost::shared_ptr<Evolver> evolver() const {return m_evolver;}
	bool canDiverge() const { return m_canDiverge; }
	PROTEIN_CHAIN_API static boost::shared_ptr<Hamiltonian> buildChain(const std::vector<double>& epsilons, double J, double kappa);
	PROTEIN_CHAIN_API static boost::shared_ptr<Hamiltonian> buildDipoles(const std::vector<double>& epsilons, double J, double kappa);
	PROTEIN_CHAIN_API static boost::shared_ptr<Hamiltonian> buildRing(const std::vector<double>& epsilons, double J, double kappa);
private:
	unsigned int m_dim;
	bool m_isHermitian;
	Eigen::MatrixXcd m_matrix;
	Eigen::VectorXcd m_eigenvalues;
	Eigen::MatrixXcd m_eigenvectors;
	boost::shared_ptr<Evolver> m_evolver;
	bool m_isNormal;
	bool m_canDiverge;
};

#endif
