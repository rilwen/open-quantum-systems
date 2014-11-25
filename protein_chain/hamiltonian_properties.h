#ifndef __HAMILTONIAN_PROPERTIES_H
#define __HAMILTONIAN_PROPERTIES_H

#include <cassert>
#include <Eigen/Core>
#include <Eigen/Eigenvalues>
#include "core.h"

//! Functions calculating properties of various hamiltonians
namespace HamiltonianProperties {
	//! Calculate eigenenergies
	template <class M> Eigen::VectorXd eigenenergies(const Eigen::MatrixXcd& hamiltonian);
	//! Calculate lowest eigenenergy
	template <class M> double lowest_eigenenergy(const Eigen::MatrixXcd& hamiltonian);
	//! Calculate highest eigenenergy
	template <class M> double highest_eigenenergy(const Eigen::MatrixXcd& hamiltonian);
	//! Convert absorption band shift C to interaction strength in the case of ring system with nearest-neighbour interactions
	PROTEIN_CHAIN_API double absorption_band_shift_to_interation_strength_ring_nn(size_t dimension, double C);
	//! Convert absorption band shift C to interaction strength in the case of wheel system with nearest-neighbour interactions
	PROTEIN_CHAIN_API double absorption_band_shift_to_interation_strength_wheel_nn(size_t dimension, double C);
};

namespace HamiltonianProperties {
	template <class M> Eigen::VectorXd eigenenergies(const M& hamiltonian)
	{
		assert(hamiltonian == hamiltonian.adjoint());
		Eigen::SelfAdjointEigenSolver<M> solver(hamiltonian);
		return solver.eigenvalues();
	}

	template<class M> double lowest_eigenenergy(const M& hamiltonian)
	{	
		return eigenenergies(hamiltonian).array().minCoeff();
	}

	template<class M> double highest_eigenenergy(const M& hamiltonian)
	{	
		return eigenenergies(hamiltonian).array().maxCoeff();
	}
}

#endif // __HAMILTONIAN_PROPERTIES_H
