#ifndef __HAMILTONIAN_FACTORY_H
#define __HAMILTONIAN_FACTORY_H

#include <Eigen/Core>
#include <vector>
#include "core.h"

//! Collection of static functions generating system Hamiltonians.
struct HamiltonianFactory
{
	PROTEIN_CHAIN_API static void buildChain(const std::vector<double>& epsilons, double J, double kappa, Eigen::MatrixXcd& hamiltonian);
	PROTEIN_CHAIN_API static void buildChain(const std::vector<double>& epsilons, double J, Eigen::MatrixXd& hamiltonian);
	PROTEIN_CHAIN_API static void buildDipoles(const std::vector<double>& epsilons, double J, double kappa, Eigen::MatrixXcd& hamiltonian);
	PROTEIN_CHAIN_API static void buildDipoles(const std::vector<double>& epsilons, double J, Eigen::MatrixXd& hamiltonian);
	PROTEIN_CHAIN_API static void buildPosDip(double epsilons, const std::vector<double>& xi, double kappa, Eigen::MatrixXcd& hamiltonian);
	PROTEIN_CHAIN_API static void buildPosDip(double epsilons, const std::vector<double>& xi, Eigen::MatrixXd& hamiltonian);
	PROTEIN_CHAIN_API static void buildWheel(const std::vector<double>& epsilons, double J, double kappa, Eigen::MatrixXcd& hamiltonian);
	PROTEIN_CHAIN_API static void buildWheel(const std::vector<double>& epsilons, double J, Eigen::MatrixXd& hamiltonian);
	PROTEIN_CHAIN_API static void buildH(const std::vector<double>& epsilons, double J, Eigen::MatrixXcd& hamiltonian);	
	PROTEIN_CHAIN_API static void buildHnobar(const std::vector<double>& epsilons, double J, Eigen::MatrixXcd& hamiltonian);	
	//! @tparam M Eigen::MatrixXd or Eigen::MatrixXcd
	template <class M> static void buildRing(const std::vector<double>& epsilons, double J, M& hamiltonian);
};

template <class M> void HamiltonianFactory::buildRing(const std::vector<double>& epsilons, double J, M& m)
{
	m.fill(0.0);

	for (unsigned int i = 0; i < epsilons.size(); ++i) {
		m(i,i) = epsilons[i];
		if (i>0) {
			m(i-1,i) = m(i,i-1) = J;
		}
	}
	if (epsilons.size() > 2)
		m(epsilons.size()-1,0) = m(0,epsilons.size()-1) = J;
}

#endif
