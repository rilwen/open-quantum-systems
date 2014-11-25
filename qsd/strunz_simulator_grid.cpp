#include "strunz_simulator_grid.h"
#include <boost/make_shared.hpp>
#include <boost/math/special_functions/fpclassify.hpp>
#include <cassert>
#include <cmath>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <vector>
#include <protein_chain/evolver_taylor_expansion.h>
#include "strunz_simulator_accumulator.h"


StrunzSimulatorGrid::StrunzSimulatorGrid(const std::vector<boost::shared_ptr<const CorrelationFunction> >& alpha, const Eigen::MatrixXcd& Hel, const double dt, const size_t nbrTimeSteps, const size_t nbrCalcStepsPerStep, const bool stochastic)
	: m_alpha(alpha), m_Hel(Hel), m_nbr_sites(alpha.size()), m_dt(dt), m_nbr_steps(nbrTimeSteps), m_nbr_calc_steps_per_step(nbrCalcStepsPerStep), m_calc_dt(dt/nbrCalcStepsPerStep), m_nbr_calc_steps(nbrTimeSteps*nbrCalcStepsPerStep)
	, m_stochastic(stochastic), m_common_alpha(false), m_Hel_norm(Hel.norm()), m_Heff(m_nbr_calc_steps + 1)
{
	if (Hel != Hel.adjoint())
		throw std::domain_error("Electronic Hamiltonian must be self-adjoint");
	if (static_cast<int>(alpha.size()) != Hel.rows())
		throw std::domain_error("Wrong alpha vector size");
	if (!nbrCalcStepsPerStep)
		throw std::domain_error("nbrCalcStepsPerStep must be at least 1");
	calculate_effective_hamiltonians_without_noise();
}

StrunzSimulatorGrid::StrunzSimulatorGrid(boost::shared_ptr<const CorrelationFunction> alpha, size_t dim, const Eigen::MatrixXcd& Hel, double dt, size_t nbrTimeSteps, const size_t nbrCalcStepsPerStep, bool stochastic)
	: m_alpha(dim, alpha), m_Hel(Hel), m_nbr_sites(dim), m_dt(dt), m_nbr_steps(nbrTimeSteps), m_nbr_calc_steps_per_step(nbrCalcStepsPerStep), m_calc_dt(dt/nbrCalcStepsPerStep), m_nbr_calc_steps(nbrTimeSteps*nbrCalcStepsPerStep)
	, m_stochastic(stochastic), m_common_alpha(true), m_Hel_norm(Hel.norm()), m_Heff(m_nbr_calc_steps + 1)
{
	if (!nbrCalcStepsPerStep)
		throw std::domain_error("nbrCalcStepsPerStep must be at least 1");
	if (Hel != Hel.adjoint())
		throw std::domain_error("Electronic Hamiltonian must be self-adjoint");
	if (static_cast<int>(dim) != Hel.rows())
		throw std::domain_error("Wrong dim");
	calculate_effective_hamiltonians_without_noise();
}

StrunzSimulatorGrid::Workspace StrunzSimulatorGrid::workspace() const
{
	return Workspace(*this);
}

double StrunzSimulatorGrid::taylor_evolver_time_step(const Eigen::MatrixXcd& effective_hamiltonian, const double base_dt) const
{
	static const size_t MAXIMUM_HAMILTONIAN_NORM_MULTIPLIER = 10;
	const double h = std::min( effective_hamiltonian.norm(), MAXIMUM_HAMILTONIAN_NORM_MULTIPLIER*m_Hel_norm );
	const double dt = std::min( h > 0 ? 0.1/h : base_dt, base_dt );
	assert(dt>0);
	return dt;
}

void StrunzSimulatorGrid::simulateTrajectory(Workspace& wksp, const Eigen::VectorXcd& initial, std::vector<Eigen::VectorXcd>& trajectory) const
{
	trajectory.resize(m_nbr_steps + 1);
	trajectory[0] = initial;	
	if (m_stochastic)
		simulate_z(wksp);
	for (size_t i = 0; i < m_nbr_steps; ++i) {
		const Eigen::MatrixXcd& heff(hamiltonian_for_evolution(wksp, i, m_dt));
		/*const double state_norm = trajectory[i].norm();
		const double avg_L = - state_norm*state_norm;*/		

		//// Euler discretization
		//trajectory[i + 1] = heff * trajectory[i];
		//trajectory[i + 1] *= std::complex<double>(0, -m_dt);
		//trajectory[i + 1] += trajectory[i];
		EvolverTaylorExpansion<4u> evolver(heff, m_dt);
		evolver.evolve(trajectory[i], m_dt, trajectory[i + 1]);
	}
}

void StrunzSimulatorGrid::calculate_effective_hamiltonians_without_noise()
{
	StrunzCalculatorGridExplicit strunz(m_alpha, m_Hel, m_nbr_calc_steps, m_calc_dt);
	StrunzCalculatorGridExplicit::Workspace wksp(strunz.workspace());
	m_Heff[0] = strunz.Hel();
	for (size_t i = 0; i < m_nbr_calc_steps; ++i) {
		strunz.step(wksp);
		m_Heff[i + 1] = strunz.effectiveHamiltonianTimesMinusI(wksp);
		m_Heff[i + 1] *= std::complex<double>(0, 1);
	}
}

void StrunzSimulatorGrid::calculate_effective_hamiltonian_with_noise(Workspace& wksp, const size_t sim_idx, double dt) const
{
	assert( dt > 0 );
	if (sim_idx < m_nbr_steps - 1) {
		wksp.m_Heff_with_noise = 0.5*m_Heff[sim_idx*m_nbr_calc_steps_per_step];
		wksp.m_Heff_with_noise += 0.5*m_Heff[(sim_idx+1)*m_nbr_calc_steps_per_step];
	} else
		wksp.m_Heff_with_noise = m_Heff[sim_idx*m_nbr_calc_steps_per_step];
	for (size_t m = 0; m < dim(); ++m) {
		if (boost::math::isnan(wksp.m_z_path(sim_idx, m)))
			throw std::runtime_error("NaN stochastic value detected");
		wksp.m_Heff_with_noise(m, m) -= wksp.m_z_path(sim_idx, m);
	}
}

const Eigen::MatrixXcd& StrunzSimulatorGrid::hamiltonian_for_evolution(Workspace& wksp, const size_t sim_idx, double dt) const
{
	if (m_stochastic) {
		calculate_effective_hamiltonian_with_noise(wksp, sim_idx, dt);
		return wksp.m_Heff_with_noise;
	} else
		return m_Heff[sim_idx*m_nbr_calc_steps_per_step];
}

void StrunzSimulatorGrid::simulate_z(Workspace& wksp) const
{
	assert(wksp.m_gauss_generator);
	assert(wksp.m_gauss_generator_wksp);
	wksp.m_gauss_generator->simulate(*wksp.m_gauss_generator_wksp, wksp.m_z_path);
}

// Nested classes

StrunzSimulatorGrid::Workspace::Workspace(const StrunzSimulatorGrid& owner)
	: m_Heff_with_noise(owner.dim(), owner.dim()), m_z_path(owner.nbrSteps(), owner.dim())
{	
	if (!owner.m_stochastic) {
		assert(!m_gauss_generator);
		assert(!m_gauss_generator_wksp);		
		owner.reset_z(*this);
	} else {
		if (owner.m_common_alpha) {
			if (owner.m_alpha.empty()) {
				boost::shared_ptr<const CorrelationFunction> empty;
				m_gauss_generator = boost::make_shared<ColoredGaussianProcessGenerator>(empty, owner.m_alpha.size(), owner.nbrSteps(), owner.dt());
			} else
				m_gauss_generator = boost::make_shared<ColoredGaussianProcessGenerator>(owner.m_alpha[0], owner.m_alpha.size(), owner.nbrSteps(), owner.dt());
		} else
			m_gauss_generator = boost::make_shared<ColoredGaussianProcessGenerator>(owner.m_alpha, owner.nbrSteps(), owner.dt());
		m_gauss_generator_wksp = boost::make_shared<ColoredGaussianProcessGenerator::Workspace>(m_gauss_generator->workspace());
	}
}
