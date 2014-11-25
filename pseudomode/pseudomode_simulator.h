#ifndef __PSEUDOMODE_SIMULATOR_H
#define __PSEUDOMODE_SIMULATOR_H

#include <vector>
#include <boost/array.hpp>
#include <boost/shared_ptr.hpp>
#include <protein_chain/evolver_sparse.h>
#include <protein_chain/reductor.h>
#include "core.h"

class ExactHamiltonian;

class PseudomodeSimulator
{
public:
	typedef boost::array<double,3> mode_param_triple;
	typedef std::vector<mode_param_triple> mode_param_row;
	typedef std::vector<mode_param_row> mode_param_array;
	PSEUDOMODE_API static const size_t LARGE_GAMMA_IDX = 0; // index of scale parameter
	PSEUDOMODE_API static const size_t SMALL_GAMMA_IDX = 1; // index of HWHM
	PSEUDOMODE_API static const size_t OMEGA_IDX = 2; // index of central frequency
	//! @param[in] modeParams Jagged array of mode params in order (Gamma_nj, gamma_nj, omega_nj)
	//! Gamma_nj is the scale parameter
	//! gamma_nj is the HWHM
	//! omega_nj is the central frequency
	PSEUDOMODE_API PseudomodeSimulator(const Eigen::MatrixXcd& Hel, const mode_param_array& modeParams, double dt, size_t nbrExcitedBathLevels);	
	PSEUDOMODE_API void convertState(const Eigen::VectorXcd& excitonState, Eigen::VectorXcd& fullState) const;
	
	//! @param[in] initial Initial state
	//! @param[in] n Number of time steps
	//! @param[out] sp Scalar products <Psi(0)|Psi^+(t)> (first element is 1 for <Psi(0)|Psi(0)>)
	PSEUDOMODE_API void simulate(const Eigen::VectorXcd& initial, size_t n, std::vector<std::complex<double> >& sp);
	
	//! @param[in] initial Initial state
	//! @param[in] n Number of time steps
	//! @param[out] rhos Reduced operators 
	//! @param[out] sp Scalar products <Psi(0)|Psi^+(t)> (first element is 1 for <Psi(0)|Psi(0)>)
	PSEUDOMODE_API void simulate(const Eigen::VectorXcd& initial, size_t n, std::vector<Eigen::MatrixXcd>& rhos, std::vector<std::complex<double> >& sp);
	
	PSEUDOMODE_API size_t dim() const { return m_dim; }
private:
	void init_omega_g(const mode_param_array& modeParams, const size_t nbr_sites, Eigen::VectorXcd& omega, Eigen::MatrixXcd& g) const;
private:	
	EvolverSparse<EigenMultiplicator> m_evolver_plus;
	boost::shared_ptr<ExactHamiltonian> m_exact_hamiltonian_generator_plus;
	EvolverSparse<EigenMultiplicator> m_evolver_minus;
	boost::shared_ptr<ExactHamiltonian> m_exact_hamiltonian_generator_minus;
	size_t m_nbr_excited_levels;
	size_t m_dim;
	Reductor m_reductor;
};

#endif // __PSEUDOMODE_SIMULATOR_H
