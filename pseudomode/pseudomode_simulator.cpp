#include "pseudomode_simulator.h"
#include <cassert>
#include <cmath>
#include <stdexcept>
#include <iostream>
#include <sstream>
#include <boost/make_shared.hpp>
#include <protein_chain/exact_hamiltonian.h>
#include <Eigen/Sparse>

void check_if_negative(double x, const char* name) {
	if (x < 0) {
		std::stringstream ss;
		ss << name << " must not be negative";
		throw std::domain_error(ss.str());
	}
}

PseudomodeSimulator::PseudomodeSimulator(const Eigen::MatrixXcd& Hel, const mode_param_array& modeParams, double dt, const size_t nbrExcitedBathLevels)
: m_evolver_plus(dt), m_evolver_minus(dt), m_nbr_excited_levels(nbrExcitedBathLevels)
{
	const size_t nbr_sites = Hel.rows();
	if (Hel.rows() != Hel.cols())
		throw std::domain_error("Hel matrix must be square");
	if (modeParams.size() != nbr_sites)
		throw std::domain_error("Wrong number of gamma rows");
	size_t mode_cnt = 0;
	for (size_t m = 0; m < nbr_sites; ++m) {
		if (modeParams[m].empty())
			throw std::domain_error("No gammas for some mode");
		mode_cnt += modeParams[m].size();
	}	
	
	Eigen::MatrixXcd g(nbr_sites, mode_cnt);	
	Eigen::VectorXcd omega(mode_cnt);
	
	init_omega_g(modeParams, nbr_sites, omega, g);	
	std::cout << g << std::endl;

	m_exact_hamiltonian_generator_plus = boost::make_shared<PredefinedElectronic>(Hel, omega, g);
	m_dim = m_exact_hamiltonian_generator_plus->dimension(nbrExcitedBathLevels);
	assert(m_dim % nbr_sites == 0);
	
	Eigen::SparseMatrix<std::complex<double> > pmHam_plus(m_dim, m_dim);
	m_exact_hamiltonian_generator_plus->generateMatrix(pmHam_plus, nbrExcitedBathLevels);
	pmHam_plus *= std::complex<double>(0, -1);
	m_evolver_plus = EvolverSparse<EigenMultiplicator>(dt, EigenMultiplicator(pmHam_plus), m_dim);
	
	m_exact_hamiltonian_generator_minus = boost::make_shared<PredefinedElectronic>(Hel, omega.conjugate(), g);
	Eigen::SparseMatrix<std::complex<double> > pmHam_minus(m_dim, m_dim);
	m_exact_hamiltonian_generator_minus->generateMatrix(pmHam_minus, nbrExcitedBathLevels);
	pmHam_minus *= std::complex<double>(0, -1);
	m_evolver_minus = EvolverSparse<EigenMultiplicator>(dt, EigenMultiplicator(pmHam_minus), m_dim);
	
	m_reductor = Reductor(nbr_sites, m_dim / nbr_sites);
}

void PseudomodeSimulator::convertState(const Eigen::VectorXcd& excitonState, Eigen::VectorXcd& fullState) const
{
	m_exact_hamiltonian_generator_plus->convertState(excitonState, fullState, m_nbr_excited_levels);
}

void PseudomodeSimulator::simulate(const Eigen::VectorXcd& initial, const size_t n, std::vector<std::complex<double> >& sp)
{
	sp.clear();
	sp.reserve(n + 1);
	sp.push_back(1.0);
	Eigen::VectorXcd state(initial);
	for (size_t i = 0; i < n; ++i) {
		m_evolver_plus.step(m_evolver_plus.timeStep(), state);
		sp.push_back(initial.dot(state));
	}
}

void PseudomodeSimulator::simulate(const Eigen::VectorXcd& initial, const size_t n, std::vector<Eigen::MatrixXcd>& rhos, std::vector<std::complex<double> >& sp)
{
	rhos.clear();
	rhos.resize(n + 1);
	sp.clear();
	sp.reserve(n + 1);
	sp.push_back(1.0);
	m_reductor.reduce(initial, rhos.front());
	Eigen::VectorXcd state_plus(initial);
	//Eigen::VectorXcd state_minus(initial);
	double cum_scale_plus = 1.0;
	for (size_t i = 0; i < n; ++i) {
		m_evolver_plus.step(m_evolver_plus.timeStep(), state_plus);
		//m_evolver_minus.step(m_evolver_minus.timeStep(), state_minus);
		const double norm_plus = state_plus.norm();
		//const double norm_minus = state_minus.norm();
		//std::cout << "NORMS: " << i*m_evolver_plus.timeStep() << " " << norm_plus << " " << norm_minus << std::endl;		
		state_plus /= norm_plus;
		//state_minus /= norm_minus;
		cum_scale_plus *= norm_plus;
		m_reductor.reduce(state_plus, state_plus, rhos[i + 1]);
		sp.push_back(initial.dot(state_plus) * cum_scale_plus);
	}
}
		
void PseudomodeSimulator::init_omega_g(const mode_param_array& modeParams, const size_t nbr_sites, Eigen::VectorXcd& omega, Eigen::MatrixXcd& g) const {
	g.setZero();
	omega.setZero();
	size_t mode_idx = 0;
	for (size_t m = 0; m < nbr_sites; ++m) {
		for (mode_param_row::const_iterator it = modeParams[m].begin(); it != modeParams[m].end(); ++it) {			
			const double Gamma = (*it)[LARGE_GAMMA_IDX];
			check_if_negative(Gamma, "Gamma");
			const double w = (*it)[OMEGA_IDX];
			check_if_negative(w, "omega");
			const double gamma = (*it)[SMALL_GAMMA_IDX];
			check_if_negative(gamma, "gamma");
			assert(mode_idx < mode_cnt);
			g(m, mode_idx) = std::complex<double>(sqrt(Gamma), 0);
			omega(mode_idx) = std::complex<double>(w, - gamma);
			++mode_idx;
		}
	}
	assert(mode_idx == mode_cnt);
}
