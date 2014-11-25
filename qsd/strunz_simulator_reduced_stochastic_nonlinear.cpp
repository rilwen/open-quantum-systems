#include "strunz_simulator_reduced_stochastic_nonlinear.h"
#include <boost/make_shared.hpp>
#include <cassert>
#include <cmath>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <vector>
#include <protein_chain/correlation_function_decomposable.h>
#include <math/ode/solver_runge_kutta_generic.h>
#include <boost/math/special_functions/fpclassify.hpp>

namespace rql {
	namespace math {
		namespace ode {
			template <> struct SolverRungeKuttaGenericNorm<StrunzSimulatorReducedStochasticNonlinear::State,true>
			{
				// norm of the (psi,integrals) state := ||psi||
				// integrated functions are proportional to psi, so if psi coverges, the integrals converge as well
				static double norm(const StrunzSimulatorReducedStochasticNonlinear::State& st) { return st.norm(); }
			};
		}		
	}
}

StrunzSimulatorReducedStochasticNonlinear::StrunzSimulatorReducedStochasticNonlinear(const std::vector<boost::shared_ptr<const CorrelationFunctionDecomposable> >& alphas, const Eigen::MatrixXcd& Hel, double dt, size_t nbr_steps)
	: m_strunz(alphas, Hel, false), m_common_alpha(false)
{
	init(dt, nbr_steps);
}

StrunzSimulatorReducedStochasticNonlinear::StrunzSimulatorReducedStochasticNonlinear(boost::shared_ptr<const CorrelationFunctionDecomposable> alpha, size_t nbr_sites, const Eigen::MatrixXcd& Hel, double dt, size_t nbr_steps)
	: m_strunz(std::vector<boost::shared_ptr<const CorrelationFunctionDecomposable> >(nbr_sites, alpha), Hel, false), m_common_alpha(true)
{
	init(dt, nbr_steps);
}

class OperatorListener
{
public:
	OperatorListener(StrunzSimulatorReducedDeterministic& strunz, std::vector<Eigen::MatrixXcd>& operators)
		: m_strunz(strunz), m_operators(operators)
	{}
	void operator()(const StrunzSimulatorReducedDeterministic::WorkspaceHamiltonian& wksp, size_t i) const
	{
		const size_t start_idx = i*m_strunz.nbr_sites();
		for (size_t m = 0; m < m_strunz.nbr_sites(); ++m) {
			// save the operators from the *end* of the i-th step
			m_strunz.calculate_operator(wksp, m, m_operators[start_idx + m]);
		}
	}
private:
	StrunzSimulatorReducedDeterministic& m_strunz;
	std::vector<Eigen::MatrixXcd>& m_operators;
};

void StrunzSimulatorReducedStochasticNonlinear::init(double dt, size_t nbr_steps)
{	
	m_dt = dt;
	m_nbr_steps = nbr_steps;
	m_nbr_hamiltonians = nbr_steps + 1;
	if (!m_strunz.nbr_sites())
		throw std::domain_error("Zero sites not supported");
	StrunzSimulatorReducedDeterministic::WorkspaceHamiltonian wksp(m_strunz.workspace_hamiltonian());
	m_hamiltonians_times_minus_i_no_noise.resize(m_nbr_hamiltonians);
	m_interaction_hamiltonians_no_noise.resize(m_nbr_hamiltonians);

	m_operators.resize(m_nbr_steps*m_strunz.nbr_sites());
	OperatorListener ubek(m_strunz, m_operators);
	m_strunz.simulate_hamiltonians_times_minus_i_adaptive(wksp, m_hamiltonians_times_minus_i_no_noise, m_dt, m_nbr_steps, ubek, &m_interaction_hamiltonians_no_noise);
	assert(m_hamiltonians_times_minus_i_no_noise.size() == m_nbr_steps+1);	
}

void StrunzSimulatorReducedStochasticNonlinear::simulate(Workspace& wksp, const Eigen::VectorXcd& initial, std::vector<Eigen::VectorXcd>& trajectory, std::vector<double>* interaction_energies) const
{
	if (initial.norm() == 0)
		throw std::domain_error("Initial state has norm zero");
	trajectory.resize(m_nbr_steps + 1);
	trajectory[0] = initial;
	wksp.reset();
	simulate_z(wksp);
	wksp.m_state.reset(initial, *this);
	Functor f(wksp);
	rql::math::ode::SolverRungeKuttaGeneric<State,false> solver(wksp.m_state);
	solver.set_norm_range(1.05);
	if (interaction_energies) {		
		// at time t=0, noise is zero in every measure
		interaction_energies->front() = initial.dot(m_interaction_hamiltonians_no_noise.front() * initial).real();
	}
	for (size_t i = 0; i < m_nbr_steps; ++i) {
		const double t = i*m_dt;
		const double next_t = (i+1)*m_dt;
		solver.solve(f, wksp.m_state, t, next_t);
		const double psi_norm = wksp.m_state.m_psi.norm();
		if (psi_norm == 0) {
			throw std::runtime_error("Zero norm");
		}
		if (!boost::math::isfinite(psi_norm)) {
			throw std::runtime_error("Norm is not finite");
		}
		wksp.m_state.m_psi /= psi_norm; // Normalisation
		trajectory[i + 1] = wksp.m_state.m_psi;
		if (interaction_energies) {
			(*interaction_energies)[i + 1] = calc_interaction_energy(wksp, next_t, wksp.m_state);
		}
	}
}

void StrunzSimulatorReducedStochasticNonlinear::simulate_z(Workspace& wksp) const
{
	assert(wksp.m_gauss_generator);
	assert(wksp.m_gauss_generator_wksp);
	wksp.m_gauss_generator->simulate(*wksp.m_gauss_generator_wksp, wksp.m_z_path);
}

StrunzSimulatorReducedStochasticNonlinear::Workspace StrunzSimulatorReducedStochasticNonlinear::workspace() const
{
	return Workspace(*this);
}

void StrunzSimulatorReducedStochasticNonlinear::calc_time_interpolation_weights_and_indices(double t, double dt, double& prev_operator_weight, size_t& prev_operator_idx, double& next_operator_weight, size_t& next_operator_idx) const
{
	prev_operator_idx = std::min(m_nbr_steps, static_cast<size_t>(floor(t/dt)));
	next_operator_idx = std::min(m_nbr_steps, static_cast<size_t>(ceil(t/dt)));
	assert(prev_operator_idx <= next_operator_idx);	
	if (prev_operator_idx == next_operator_idx) {
		prev_operator_weight = 1;
		next_operator_weight = 0;
	} else {
		assert(next_operator_idx == prev_operator_idx + 1);
		prev_operator_weight = (t - prev_operator_idx*dt)/dt;
		next_operator_weight = (next_operator_idx*dt - t)/dt;
	}
}

StrunzSimulatorReducedStochasticNonlinear::InterpolationParams StrunzSimulatorReducedStochasticNonlinear::interpolate_calculator_outputs(Workspace& wksp, double t, const std::vector<Eigen::MatrixXcd>& outputs, Eigen::MatrixXcd& result) const
{
	InterpolationParams ip;
	calc_time_interpolation_weights_and_indices(t, m_dt, ip.weight_prev, ip.prev_ham_idx, ip.weight_next, ip.next_ham_idx);
	ip.noise_idx = std::min(wksp.m_owner.m_nbr_steps - 1, ip.prev_ham_idx);
	ip.prev_operator_start_idx = (ip.prev_ham_idx - 1)*m_strunz.nbr_sites();
	ip.next_operator_start_idx = (ip.next_ham_idx - 1)*m_strunz.nbr_sites();

	result = ip.weight_prev*outputs[ip.prev_ham_idx];
	if (ip.weight_next) {
		result += ip.weight_next*outputs[ip.next_ham_idx];
	}
	return ip;
}

double StrunzSimulatorReducedStochasticNonlinear::calc_interaction_energy(Workspace& wksp, double t, const State& y) const
{
	const double y_norm = y.m_psi.norm();
	const double y_norm_sqr = y_norm*y_norm;
	if (y_norm < 0.9) {
		std::cerr << "Norm is not conserved: " << y_norm << " at time " << t << std::endl;
	}
	const rql::math::Jagged2DArray<std::complex<double> >& y_drift_integrals = y.m_drift_integrals;

	Eigen::MatrixXcd& ih = wksp.m_interaction_hamiltonian;
	const InterpolationParams interpolation_params = interpolate_calculator_outputs(wksp, t, m_interaction_hamiltonians_no_noise, ih);
	/*size_t prev_ham_idx;
	size_t next_ham_idx;
	double weight_prev;
	double weight_next; 
	wksp.m_owner.calc_time_interpolation_weights_and_indices(t, m_dt, weight_prev, prev_ham_idx, weight_next, next_ham_idx);
	const size_t noise_idx = std::min(wksp.m_owner.m_nbr_steps - 1, prev_ham_idx);
	const size_t prev_operator_start_idx = (prev_ham_idx - 1)*m_strunz.nbr_sites();
	const size_t next_operator_start_idx = (next_ham_idx - 1)*m_strunz.nbr_sites();

	ih = weight_prev*m_interaction_hamiltonians_no_noise[prev_ham_idx];
	if (weight_next) {
		ih += weight_next*m_interaction_hamiltonians_no_noise[next_ham_idx];
	}*/
	
	Eigen::VectorXcd& tmp_vec = wksp.m_tmp_vec;
	const std::vector<Eigen::MatrixXcd>& operators = wksp.m_owner.m_operators;
	const Eigen::VectorXcd& y_psi = y.m_psi;
	for (size_t m = 0; m < m_strunz.nbr_sites(); ++m) {		
		const std::complex<double>& y_psi_m = y_psi[m];
		const double avg_L = - (y_psi_m.real()*y_psi_m.real() + y_psi_m.imag()*y_psi_m.imag()) / y_norm_sqr; // <L_m>_t
		const size_t nbr_peaks = y.m_drift_integrals.row_size(m);
		std::complex<double> noise_with_drift = wksp.m_z_path(interpolation_params.noise_idx, m);
		if (!boost::math::isfinite(noise_with_drift.real()) || !boost::math::isfinite(noise_with_drift.imag())) {
			throw std::runtime_error("Noise is not finite");
		}
		for (size_t k = 0; k < nbr_peaks; ++k) {
			const std::complex<double>& ymk = y_drift_integrals(m, k);
			noise_with_drift += ymk;
		}
		ih(m,m) -= noise_with_drift * std::complex<double>(0, 1);
		// add the term - <L_m^dagger>_t D_m(t)
		// D_m(t) is calculated by interpolation
		if (interpolation_params.prev_ham_idx > 0) { // at time 0, all D_m == 0			
			ih -= (avg_L * interpolation_params.weight_prev * std::complex<double>(0, 1)) * operators[interpolation_params.prev_operator_start_idx + m];
		}
		if (interpolation_params.weight_next && interpolation_params.next_ham_idx > 0) {
			ih -= (avg_L * interpolation_params.weight_next * std::complex<double>(0, 1)) * operators[interpolation_params.next_operator_start_idx + m];
		}
	}

	return y.m_psi.dot(ih * y.m_psi).real();
}

// Nested classes

StrunzSimulatorReducedStochasticNonlinear::State::State(const StrunzSimulatorReducedStochasticNonlinear& owner)
	: m_psi(owner.m_strunz.nbr_sites())
{
	std::vector<size_t> pc(owner.m_strunz.correlation_functions_decomposition().peak_counts());
	m_drift_integrals = rql::math::Jagged2DArray<std::complex<double> >(pc.begin(), pc.end());
}

void StrunzSimulatorReducedStochasticNonlinear::State::reset(const Eigen::VectorXcd& initial, const StrunzSimulatorReducedStochasticNonlinear& owner)
{
	m_psi = initial;	
	std::fill(m_drift_integrals.flat_begin(), m_drift_integrals.flat_end(), 0.0);
}

StrunzSimulatorReducedStochasticNonlinear::State& StrunzSimulatorReducedStochasticNonlinear::State::operator+=(const State& rhs)
{
	m_psi += rhs.m_psi;
	const rql::math::Jagged2DArray<std::complex<double> >::const_flat_iterator l_end = m_drift_integrals.flat_end();
	rql::math::Jagged2DArray<std::complex<double> >::const_flat_iterator r_it = rhs.m_drift_integrals.flat_begin();
	for (rql::math::Jagged2DArray<std::complex<double> >::flat_iterator l_it = m_drift_integrals.flat_begin(); l_it != l_end; ++l_it,++r_it) {
		assert(r_it != rhs.m_drift_integrals.flat_end());
		(*l_it) += *r_it;
	}
	return *this;
}

StrunzSimulatorReducedStochasticNonlinear::State& StrunzSimulatorReducedStochasticNonlinear::State::operator-=(const State& rhs)
{
	m_psi -= rhs.m_psi;
	const rql::math::Jagged2DArray<std::complex<double> >::const_flat_iterator l_end = m_drift_integrals.flat_end();
	rql::math::Jagged2DArray<std::complex<double> >::const_flat_iterator r_it = rhs.m_drift_integrals.flat_begin();
	for (rql::math::Jagged2DArray<std::complex<double> >::flat_iterator l_it = m_drift_integrals.flat_begin(); l_it != l_end; ++l_it,++r_it) {
		assert(r_it != rhs.m_drift_integrals.flat_end());
		(*l_it) -= *r_it;
	}
	return *this;
}

StrunzSimulatorReducedStochasticNonlinear::State& StrunzSimulatorReducedStochasticNonlinear::State::operator*=(const std::complex<double>& x)
{
	m_psi *= x;
	const rql::math::Jagged2DArray<std::complex<double> >::const_flat_iterator l_end = m_drift_integrals.flat_end();
	for (rql::math::Jagged2DArray<std::complex<double> >::flat_iterator l_it = m_drift_integrals.flat_begin(); l_it != l_end; ++l_it) {
		(*l_it) *= x;
	}
	return *this;
}

StrunzSimulatorReducedStochasticNonlinear::State& StrunzSimulatorReducedStochasticNonlinear::State::operator/=(const std::complex<double>& x)
{
	m_psi /= x;
	const rql::math::Jagged2DArray<std::complex<double> >::const_flat_iterator l_end = m_drift_integrals.flat_end();
	for (rql::math::Jagged2DArray<std::complex<double> >::flat_iterator l_it = m_drift_integrals.flat_begin(); l_it != l_end; ++l_it) {
		(*l_it) /= x;
	}
	return *this;
}

StrunzSimulatorReducedStochasticNonlinear::Workspace::Workspace(const StrunzSimulatorReducedStochasticNonlinear& owner)
	: m_owner(owner), m_z_path(owner.m_nbr_steps, owner.m_strunz.nbr_sites())
	, m_nonlinear_hamiltonian_times_i(owner.m_strunz.nbr_sites(), owner.m_strunz.nbr_sites())
	, m_nonlinear_hamiltonian(owner.m_strunz.nbr_sites(), owner.m_strunz.nbr_sites())
	, m_interaction_hamiltonian(owner.m_strunz.nbr_sites(), owner.m_strunz.nbr_sites())
	, m_state(owner)
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

void StrunzSimulatorReducedStochasticNonlinear::Workspace::reset()
{
	m_z_path.setZero();
}

StrunzSimulatorReducedStochasticNonlinear::Functor::Functor(Workspace& wksp)
	: m_nbr_sites(wksp.m_owner.m_strunz.nbr_sites()), m_wksp(wksp), m_alpha_decomp(wksp.m_owner.m_strunz.correlation_functions_decomposition())
{
}

void StrunzSimulatorReducedStochasticNonlinear::Functor::operator()(double t, const State& y, State& dy)
{
	const double y_norm = y.m_psi.norm();
	const double y_norm_sqr = y_norm*y_norm;
	if (y_norm < 0.9) {
		std::cerr << "Norm is not conserved: " << y_norm << " at time " << t << std::endl;
	}
	const rql::math::Jagged2DArray<std::complex<double> >& y_drift_integrals = y.m_drift_integrals;
	rql::math::Jagged2DArray<std::complex<double> >& dy_drift_integrals = dy.m_drift_integrals;

	const double dt = m_wksp.m_owner.m_dt;
	assert(dt>0);
	Eigen::MatrixXcd& hti = m_wksp.m_nonlinear_hamiltonian_times_i;
	const InterpolationParams interpolation_params = m_wksp.m_owner.interpolate_calculator_outputs(m_wksp, t, m_wksp.m_owner.m_hamiltonians_times_minus_i_no_noise, hti);
	/*size_t prev_ham_idx;
	size_t next_ham_idx;
	double weight_prev;
	double weight_next; 
	m_wksp.m_owner.calc_time_interpolation_weights_and_indices(t, dt, weight_prev, prev_ham_idx, weight_next, next_ham_idx);
	const size_t noise_idx = std::min(m_wksp.m_owner.m_nbr_steps - 1, prev_ham_idx);
	const size_t prev_operator_start_idx = (prev_ham_idx - 1)*m_nbr_sites;
	const size_t next_operator_start_idx = (next_ham_idx - 1)*m_nbr_sites;

	hti = weight_prev*m_wksp.m_owner.m_hamiltonians_times_minus_i_no_noise[prev_ham_idx];
	if (weight_next) {
		hti += weight_next*m_wksp.m_owner.m_hamiltonians_times_minus_i_no_noise[next_ham_idx];
	}*/
	
	Eigen::VectorXcd& tmp_vec = m_wksp.m_tmp_vec;
	const std::vector<Eigen::MatrixXcd>& operators = m_wksp.m_owner.m_operators;
	const Eigen::VectorXcd& y_psi = y.m_psi;
	for (size_t m = 0; m < m_nbr_sites; ++m) {		
		const std::complex<double>& y_psi_m = y_psi[m];
		const double avg_L = - (y_psi_m.real()*y_psi_m.real() + y_psi_m.imag()*y_psi_m.imag()) / y_norm_sqr; // <L_m>_t
		const size_t nbr_peaks = y.m_drift_integrals.row_size(m);
		std::complex<double> noise_with_drift = m_wksp.m_z_path(interpolation_params.noise_idx, m);
		if (!boost::math::isfinite(noise_with_drift.real()) || !boost::math::isfinite(noise_with_drift.imag())) {
			throw std::runtime_error("Noise is not finite");
		}
		for (size_t k = 0; k < nbr_peaks; ++k) {
			const std::complex<double>& ymk = y_drift_integrals(m, k);
			dy_drift_integrals(m,k) = m_alpha_decomp.scale(m, k)*avg_L + conj(m_alpha_decomp.exponent(m, k))*ymk;
			noise_with_drift += ymk;
		}
		hti(m,m) -= noise_with_drift;
		// add the term - <L_m^dagger>_t D_m(t)
		// D_m(t) is calculated by interpolation
		if (interpolation_params.prev_ham_idx > 0) { // at time 0, all D_m == 0			
			hti -= (avg_L * interpolation_params.weight_prev) * operators[interpolation_params.prev_operator_start_idx + m];
		}
		if (interpolation_params.weight_next && interpolation_params.next_ham_idx > 0) {
			hti -= (avg_L * interpolation_params.weight_next) * operators[interpolation_params.next_operator_start_idx + m];
		}
	}

	std::complex<double> a = y_psi.dot(hti*y_psi);
	a += y_psi.dot(hti.adjoint()*y_psi);
	a /= 2*y_norm_sqr;

	dy.m_psi = hti*y_psi;
	dy.m_psi -= a*y_psi;
}

