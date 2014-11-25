#include "strunz_simulator_reduced_stochastic_linear.h"
#include <boost/make_shared.hpp>
#include <cassert>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <vector>
#include <protein_chain/correlation_function_decomposable.h>


StrunzSimulatorReducedStochasticLinear::StrunzSimulatorReducedStochasticLinear(const std::vector<boost::shared_ptr<const CorrelationFunctionDecomposable> >& alphas, const Eigen::MatrixXcd& Hel, double dt, size_t nbr_steps)
	: m_strunz(alphas, Hel, true), m_dt(dt), m_nbr_steps(nbr_steps), m_common_alpha(false)
{
	init();
}

StrunzSimulatorReducedStochasticLinear::StrunzSimulatorReducedStochasticLinear(boost::shared_ptr<const CorrelationFunctionDecomposable> alpha, size_t nbr_sites, const Eigen::MatrixXcd& Hel, double dt, size_t nbr_steps)
	: m_strunz(std::vector<boost::shared_ptr<const CorrelationFunctionDecomposable> >(nbr_sites, alpha), Hel, true), m_dt(dt), m_nbr_steps(nbr_steps), m_common_alpha(true)
{
	init();
}

void StrunzSimulatorReducedStochasticLinear::init()
{	
	if (!m_strunz.nbr_sites())
		throw std::domain_error("Zero sites not supported");
	StrunzSimulatorReducedDeterministic::WorkspaceHamiltonian wksp(m_strunz.workspace_hamiltonian());
	m_hamiltonians_times_minus_i_no_noise.resize(4 * m_nbr_steps);
	std::vector<Eigen::MatrixXcd> hamiltonians;

	m_strunz.simulate_hamiltonians_times_minus_i_adaptive(wksp, hamiltonians, m_dt/2, 2*m_nbr_steps);

	assert(hamiltonians.size() == 2*m_nbr_steps+1);	
	for (size_t i = 0; i < m_nbr_steps; ++i) {
		m_hamiltonians_times_minus_i_no_noise[4*i] = hamiltonians[2*i];
		m_hamiltonians_times_minus_i_no_noise[4*i+1] = hamiltonians[2*i+1];
		m_hamiltonians_times_minus_i_no_noise[4*i+2] = hamiltonians[2*i+1];
		m_hamiltonians_times_minus_i_no_noise[4*i+3] = hamiltonians[2*i+2];
	}
}

void StrunzSimulatorReducedStochasticLinear::simulate(Workspace& wksp, const Eigen::VectorXcd& initial, std::vector<Eigen::VectorXcd>& trajectory) const
{
	trajectory.resize(m_nbr_steps + 1);
	trajectory[0] = initial;
	wksp.reset();
	simulate_z(wksp);
	calculate_hamiltonians_with_noise(wksp);
	m_strunz.simulate_psi(wksp.m_strunz_wksp_state, wksp.m_hamiltonians_times_minus_i_with_noise, initial, trajectory, m_dt, m_nbr_steps);
}

void StrunzSimulatorReducedStochasticLinear::simulate_z(Workspace& wksp) const
{
	assert(wksp.m_gauss_generator);
	assert(wksp.m_gauss_generator_wksp);
	wksp.m_gauss_generator->simulate(*wksp.m_gauss_generator_wksp, wksp.m_z_path);
}

void StrunzSimulatorReducedStochasticLinear::calculate_hamiltonians_with_noise(Workspace& wksp) const
{
	assert( m_hamiltonians_times_minus_i_no_noise.size() == wksp.m_hamiltonians_times_minus_i_with_noise.size() );
	std::copy(m_hamiltonians_times_minus_i_no_noise.begin(), m_hamiltonians_times_minus_i_no_noise.end(), wksp.m_hamiltonians_times_minus_i_with_noise.begin());
	for (size_t i = 0; i < m_nbr_steps; ++i) {
		const size_t j = 4*i;
		for (size_t m = 0; m < m_strunz.nbr_sites(); ++m) {
			for (size_t k = 0; k < 4; ++k)
				wksp.m_hamiltonians_times_minus_i_with_noise[j + k](m, m) -= wksp.m_z_path(i, m);
		}
	}
}

StrunzSimulatorReducedStochasticLinear::Workspace StrunzSimulatorReducedStochasticLinear::workspace() const
{
	return Workspace(*this);
}

// Nested classes

StrunzSimulatorReducedStochasticLinear::Workspace::Workspace(const StrunzSimulatorReducedStochasticLinear& owner)
	: m_z_path(owner.m_nbr_steps, owner.m_strunz.nbr_sites()), m_strunz_wksp_state(owner.m_strunz.workspace_state())
	, m_hamiltonians_times_minus_i_with_noise(owner.m_hamiltonians_times_minus_i_no_noise.size())
{	
	assert(!owner.m_strunz.alphas().empty());
	std::vector<boost::shared_ptr<const CorrelationFunction> > alphas_abstract(owner.m_strunz.alphas().size());
	std::copy(owner.m_strunz.alphas().begin(), owner.m_strunz.alphas().end(), alphas_abstract.begin());
	
	if (owner.m_common_alpha) {
		m_gauss_generator = boost::make_shared<ColoredGaussianProcessGenerator>(alphas_abstract[0], alphas_abstract.size(), owner.m_nbr_steps, owner.m_dt);
	} else
		m_gauss_generator = boost::make_shared<ColoredGaussianProcessGenerator>(alphas_abstract, owner.m_nbr_steps, owner.m_dt);
	m_gauss_generator_wksp = boost::make_shared<ColoredGaussianProcessGenerator::Workspace>(m_gauss_generator->workspace());
}

void StrunzSimulatorReducedStochasticLinear::Workspace::reset()
{
	m_z_path.setZero();
}




