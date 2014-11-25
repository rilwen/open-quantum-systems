#include "strunz_simulator_reduced_deterministic.h"
#include "correlation_functions_decomposition_factory.h"
#include <cassert>
#include <cmath>
#include <iostream>
#include <stdexcept>

StrunzSimulatorReducedDeterministic::StrunzSimulatorReducedDeterministic(const std::vector<boost::shared_ptr<const CorrelationFunctionDecomposable> >& alphas, const Eigen::MatrixXcd& Hel, bool correct_psi_norm)
	: m_calc(CorrelationFunctionsDecompositionFactory::decomposeVirtual(alphas), Hel), m_alphas(alphas), m_correct_psi_norm(correct_psi_norm)
{
	if (alphas.size() != Hel.rows() || Hel.rows() != Hel.cols())
		throw std::domain_error("Bad inputs");
}

void StrunzSimulatorReducedDeterministic::simulate(Workspace& wksp, const Eigen::VectorXcd& initial, std::vector<Eigen::VectorXcd>& trajectory, const double dt, const size_t nbr_steps, Diagnostics* diagnostics) const
{
	if (initial.size() != m_calc.nbrSites())
		throw std::domain_error("Bad input dimension");
	if (dt <= 0)
		throw std::domain_error("Dt must be positive");
	trajectory.resize(nbr_steps + 1);
	trajectory[0] = initial;
	wksp.clear();
	for (size_t i = 0; i < nbr_steps; ++i) {
		const double t = i*dt;
		step_both(wksp, trajectory[i], trajectory[i + 1], t, dt);
		if (diagnostics) {
			diagnostics->add_trace_deviation(t + dt, m_calc.total_absolute_trace_deviation_from_theoretical(wksp.hamiltonian_workspace().m_operators, t+dt));
			diagnostics->add_hamiltonian_norm(t, wksp.hamiltonian_workspace().m_h1.norm());
		}
	}
}

void StrunzSimulatorReducedDeterministic::simulate(Workspace& wksp, const Eigen::VectorXcd& initial, std::vector<std::complex<double> >& scalar_products, const double dt, const size_t nbr_steps) const
{
	if (initial.size() != m_calc.nbrSites())
		throw std::domain_error("Bad input dimension");
	if (dt <= 0)
		throw std::domain_error("Dt must be positive");
	scalar_products.resize(nbr_steps + 1);
	scalar_products[0] = initial.dot(initial);
	wksp.clear();
	Eigen::VectorXcd prev = initial;
	Eigen::VectorXcd next;
	for (size_t i = 0; i < nbr_steps; ++i) {
		step_both(wksp, prev, next, i*dt, dt);
		scalar_products[i + 1] = initial.dot(next);
		prev = next;
	}
}

template <class O> static void update_min_max_norm(const O& algebraic_object, double& min_norm, double& max_norm)
{
	min_norm = std::min(algebraic_object.norm(), min_norm);
	max_norm = std::max(algebraic_object.norm(), max_norm);
}

bool StrunzSimulatorReducedDeterministic::norms_diverged(double min_norm, double max_norm)
{
	assert(min_norm >= 0);
	assert(max_norm >= 0);
	return min_norm > 0 && max_norm >= MAXIMUM_MAX_MIN_NORM_RATIO*min_norm;
}

bool StrunzSimulatorReducedDeterministic::step_hamiltonian_impl(WorkspaceHamiltonian& wksp, const double t, const double dt, const bool force) const
{
	const StrunzCalculatorReduced::Data& prev_state = wksp.m_operators;	

	// 4th order Runge-Kutta for operator evolution
	m_calc.calculate_effective_hamiltonian_times_minus_i_and_time_derivative(prev_state, wksp.m_h1, wksp.m_k1);
	wksp.m_k1 *= dt;
	double min_h_norm = wksp.m_h1.norm();
	double max_h_norm = min_h_norm;

	wksp.m_operators_next.copy_data_from(prev_state);
	wksp.m_operators_next.add(wksp.m_k1, 0.5);
	//m_calc.normalize_traces(wksp.m_operators_next, t + 0.5*dt);
	m_calc.calculate_effective_hamiltonian_times_minus_i_and_time_derivative(wksp.m_operators_next, wksp.m_h2, wksp.m_k2);
	wksp.m_k2 *= dt;
	update_min_max_norm(wksp.m_h2, min_h_norm, max_h_norm);
	if (norms_diverged(min_h_norm, max_h_norm) && !force)
		return false;

	wksp.m_operators_next.copy_data_from(prev_state);
	wksp.m_operators_next.add(wksp.m_k2, 0.5);
	//m_calc.normalize_traces(wksp.m_operators_next, t + 0.5*dt);
	m_calc.calculate_effective_hamiltonian_times_minus_i_and_time_derivative(wksp.m_operators_next, wksp.m_h3, wksp.m_k3);
	wksp.m_k3 *= dt;
	update_min_max_norm(wksp.m_h3, min_h_norm, max_h_norm);
	if (norms_diverged(min_h_norm, max_h_norm) && !force)
		return false;

	wksp.m_operators_next.copy_data_from(prev_state);
	wksp.m_operators_next += wksp.m_k3;
	//m_calc.normalize_traces(wksp.m_operators_next, t + dt);
	m_calc.calculate_effective_hamiltonian_times_minus_i_and_time_derivative(wksp.m_operators_next, wksp.m_h4, wksp.m_k4);
	wksp.m_k4 *= dt;
	update_min_max_norm(wksp.m_h4, min_h_norm, max_h_norm);
	if (norms_diverged(min_h_norm, max_h_norm) && !force)
		return false;

	wksp.m_operators_next.copy_data_from(wksp.m_k2);
	wksp.m_operators_next += wksp.m_k3;
	wksp.m_operators_next *= 2.0;
	wksp.m_operators_next += wksp.m_k1;
	wksp.m_operators_next += wksp.m_k4;
	wksp.m_operators_next /= 6.0;
	wksp.m_operators += wksp.m_operators_next;
	//m_calc.normalize_traces(wksp.m_operators, t + dt);
	return true;
}

void StrunzSimulatorReducedDeterministic::step_hamiltonian_adaptive(WorkspaceHamiltonian& wksp, const double t, const double dt, const unsigned int nbr_allowed_dt_divisions) const
{
	if (nbr_allowed_dt_divisions == 0) {
		step_hamiltonian_impl(wksp, t, dt, true);
	} else {
		if (!step_hamiltonian_impl(wksp, t, dt, false)) {
			assert( nbr_allowed_dt_divisions > 0 );
			const double new_dt = dt / 2;
			const unsigned int new_nbr_allowed_dt_divisions = nbr_allowed_dt_divisions - 1;
			step_hamiltonian_adaptive(wksp, t, new_dt, new_nbr_allowed_dt_divisions);
			step_hamiltonian_adaptive(wksp, t + new_dt, new_dt, new_nbr_allowed_dt_divisions);
		}
	}
}

static void print_details(const char* name, const Eigen::MatrixXcd& h_times_minus_i)
{
	std::cerr << name << std::endl;
	std::cerr << "Norm\tNormRe\tNormIm\n";
	std::cerr << h_times_minus_i.norm() << "\t" << 0.5*(h_times_minus_i + h_times_minus_i.adjoint()).norm() << "\t" << 0.5*(h_times_minus_i - h_times_minus_i.adjoint()).norm() << std::endl;
}

// delta and next can be the same vectors
static void psi_add_delta_with_correction(const Eigen::VectorXcd& prev, const Eigen::VectorXcd& delta, const double log_norm_squared_delta, const double scale_factor, Eigen::VectorXcd& next, bool apply_correction)
{
	next = delta;
	if (scale_factor != 1)
		next *= scale_factor;
	next += prev;
	if (apply_correction) {
		// calculate correct norm for next psi (assumes that prev norm was correct)
		const double old_norm_squared = next.dot(next).real();
		if (old_norm_squared > 0) {
			const double new_norm_squared = (prev.dot(prev).real()) * exp(scale_factor*log_norm_squared_delta);
			assert(new_norm_squared >= 0);
			const double correction_factor = sqrt(new_norm_squared / old_norm_squared);
			next *= correction_factor; // comment out this line to disable corrections
		}
	}
}

bool StrunzSimulatorReducedDeterministic::step_psi_impl(const Eigen::MatrixXcd& h1, const Eigen::MatrixXcd& h2, const Eigen::MatrixXcd& h3, const Eigen::MatrixXcd& h4, WorkspaceState& wksp_st, const Eigen::VectorXcd& prev, Eigen::VectorXcd& next, const double dt, const bool force) const
{		
	// 4th order Runge-Kutta for state evolution
	wksp_st.m_psi_k1.noalias() = h1*prev;		
	wksp_st.m_psi_k1 *= dt;
	double min_k_norm = wksp_st.m_psi_k1.norm();
	double max_k_norm = min_k_norm;
	const double lognorm_k1 = dt * calculate_log_norm_derivative(wksp_st, prev, h1);

	psi_add_delta_with_correction(prev, wksp_st.m_psi_k1, lognorm_k1, 0.5, next, m_correct_psi_norm);
	wksp_st.m_psi_k2.noalias() = h2*next;
	wksp_st.m_psi_k2 *= dt;
	update_min_max_norm(wksp_st.m_psi_k2, min_k_norm, max_k_norm);
	const double lognorm_k2 = dt * calculate_log_norm_derivative(wksp_st, next, h2);
	if (norms_diverged(min_k_norm, max_k_norm) && !force)
		return false;

	psi_add_delta_with_correction(prev, wksp_st.m_psi_k2, lognorm_k2, 0.5, next, m_correct_psi_norm);
	wksp_st.m_psi_k3.noalias() = h3*next;
	wksp_st.m_psi_k3 *= dt;
	update_min_max_norm(wksp_st.m_psi_k3, min_k_norm, max_k_norm);
	if (norms_diverged(min_k_norm, max_k_norm) && !force)
		return false;
	const double lognorm_k3 = dt * calculate_log_norm_derivative(wksp_st, next, h3);

	psi_add_delta_with_correction(prev, wksp_st.m_psi_k3, lognorm_k3, 1.0, next, m_correct_psi_norm);
	wksp_st.m_psi_k4.noalias() = h4*next;
	wksp_st.m_psi_k4 *= dt;
	update_min_max_norm(wksp_st.m_psi_k4, min_k_norm, max_k_norm);
	if (norms_diverged(min_k_norm, max_k_norm) && !force)
		return false;
	const double lognorm_k4 = dt * calculate_log_norm_derivative(wksp_st, next, h4);

	next = wksp_st.m_psi_k2;
	next += wksp_st.m_psi_k3;
	next *= 2;
	next += wksp_st.m_psi_k1;
	next += wksp_st.m_psi_k4;
	next /= 6;
	psi_add_delta_with_correction(prev, next, (lognorm_k1 + 2*(lognorm_k2 + lognorm_k3) + lognorm_k4)/6.0, 1.0, next, m_correct_psi_norm);	
	return true;
}

void StrunzSimulatorReducedDeterministic::step_both_adaptive(WorkspaceHamiltonian& wksp_ham, WorkspaceState& wksp_st, const Eigen::VectorXcd& prev, Eigen::VectorXcd& next, const double t, const double dt, const unsigned int nbr_allowed_dt_divisions) const
{	
	if (nbr_allowed_dt_divisions == 0) {
		// just do it
		step_hamiltonian_impl(wksp_ham, t, dt, true);
		step_psi_adaptive(wksp_ham, wksp_st, prev, next, dt);
	} else {
		const bool step_ok = step_hamiltonian_impl(wksp_ham, t, dt, false);
		if (!step_ok) {
			assert( nbr_allowed_dt_divisions > 0 );
			// take two shorter steps
			const double new_dt = dt / 2;
			const unsigned int new_nbr_allowed_dt_divisions = nbr_allowed_dt_divisions - 1;
			step_both_adaptive(wksp_ham, wksp_st, prev, wksp_st.m_mid_psi[new_nbr_allowed_dt_divisions], t, new_dt, new_nbr_allowed_dt_divisions);
			step_both_adaptive(wksp_ham, wksp_st, wksp_st.m_mid_psi[new_nbr_allowed_dt_divisions], next, t + new_dt, new_dt, new_nbr_allowed_dt_divisions);
		} else {
			step_psi_adaptive(wksp_ham, wksp_st, prev, next, dt);
		}
	}
}

void StrunzSimulatorReducedDeterministic::step_psi_adaptive(const Eigen::MatrixXcd& h1, const Eigen::MatrixXcd& h2, const Eigen::MatrixXcd& h3, const Eigen::MatrixXcd& h4, WorkspaceState& wksp_st, const Eigen::VectorXcd& prev, Eigen::VectorXcd& next, const double dt, const unsigned int nbr_allowed_dt_divisions) const
{
	if (nbr_allowed_dt_divisions == 0) {
		// just do it
		step_psi_impl(h1, h2, h3, h4, wksp_st, prev, next, dt, true);
	} else {
		const bool step_ok = step_psi_impl(h1, h2, h3, h4, wksp_st, prev, next, dt, false);
		if (!step_ok) {
			assert(nbr_allowed_dt_divisions>0);
			// take two shorter steps
			const double new_dt = dt / 2;
			const unsigned int new_nbr_allowed_dt_divisions = nbr_allowed_dt_divisions - 1;
			wksp_st.m_mid_h[new_nbr_allowed_dt_divisions] = 0.5*h1;
			wksp_st.m_mid_h[new_nbr_allowed_dt_divisions] += 0.5*h2;
			step_psi_adaptive(h1, wksp_st.m_mid_h[new_nbr_allowed_dt_divisions], wksp_st.m_mid_h[new_nbr_allowed_dt_divisions], h2, wksp_st, prev, wksp_st.m_mid_psi[new_nbr_allowed_dt_divisions], new_dt, new_nbr_allowed_dt_divisions);
			wksp_st.m_mid_h[new_nbr_allowed_dt_divisions] = 0.5*h3;
			wksp_st.m_mid_h[new_nbr_allowed_dt_divisions] += 0.5*h4;
			step_psi_adaptive(h3, wksp_st.m_mid_h[new_nbr_allowed_dt_divisions], wksp_st.m_mid_h[new_nbr_allowed_dt_divisions], h4, wksp_st, wksp_st.m_mid_psi[new_nbr_allowed_dt_divisions], next, new_dt, new_nbr_allowed_dt_divisions);
		}
	}
}

void StrunzSimulatorReducedDeterministic::simulate_hamiltonians_times_minus_i_nonadaptive(WorkspaceHamiltonian& wksp, std::vector<Eigen::MatrixXcd>& hamiltonians_times_minus_i, const double dt, const size_t nbr_steps, const bool all) const
{
	wksp.clear();
	if (all)
		hamiltonians_times_minus_i.resize(4*nbr_steps + 1);
	else
		hamiltonians_times_minus_i.resize(nbr_steps + 1);
	for (size_t i = 0; i < nbr_steps; ++i) {
		step_hamiltonian_impl(wksp, i*dt, dt, true); // force a single RK4 step
		if (all) {
			hamiltonians_times_minus_i[4*i] = wksp.m_h1;
			hamiltonians_times_minus_i[4*i + 1] = wksp.m_h2;
			hamiltonians_times_minus_i[4*i + 2] = wksp.m_h3;
			hamiltonians_times_minus_i[4*i + 3] = wksp.m_h4;
		} else {
			hamiltonians_times_minus_i[i] = wksp.m_h1;
		}
	}
	// save the last one
	m_calc.calculate_effective_hamiltonian_times_minus_i(wksp.m_operators, hamiltonians_times_minus_i.back());
}

struct ListenerNoOp
{
	void operator()(const StrunzSimulatorReducedDeterministic::WorkspaceHamiltonian& wksp, size_t i) const {}
};

void StrunzSimulatorReducedDeterministic::simulate_hamiltonians_times_minus_i_adaptive(WorkspaceHamiltonian& wksp, std::vector<Eigen::MatrixXcd>& hamiltonians_times_minus_i, const double dt, const size_t nbr_steps, std::vector<Eigen::MatrixXcd>* interaction_hamiltonians) const
{
	ListenerNoOp list_no_op;
	simulate_hamiltonians_times_minus_i_adaptive(wksp, hamiltonians_times_minus_i, dt, nbr_steps, list_no_op, interaction_hamiltonians);
}

void StrunzSimulatorReducedDeterministic::simulate_psi(WorkspaceState& wksp, const std::vector<Eigen::MatrixXcd>& hamiltonians_times_minus_i, const Eigen::VectorXcd& initial, std::vector<Eigen::VectorXcd>& trajectory, double dt, const size_t nbr_steps) const
{
	if (hamiltonians_times_minus_i.size() < 4*nbr_steps)
		throw std::domain_error("Not enough hamiltonians provided");
	trajectory.resize(nbr_steps + 1);
	trajectory[0] = initial;
	for (size_t i = 0; i < nbr_steps; ++i) {
		step_psi_adaptive(hamiltonians_times_minus_i[4*i], hamiltonians_times_minus_i[4*i+1], hamiltonians_times_minus_i[4*i+2], hamiltonians_times_minus_i[4*i+3], wksp, trajectory[i], trajectory[i + 1], dt, MAX_ALLOWED_DT_DIVISIONS);
	}
}

double StrunzSimulatorReducedDeterministic::calculate_log_norm_derivative(const Eigen::VectorXcd& psi, const Eigen::MatrixXcd& h_eff_times_minus_i, Eigen::MatrixXcd& tmp) const
{
	const double norm_squared = psi.dot(psi).real();
	assert(norm_squared >= 0);
	if (norm_squared == 0)
		return 0;
	tmp = h_eff_times_minus_i;
	tmp += h_eff_times_minus_i.adjoint();
	assert( tmp == tmp.adjoint() );
	const double mean_value = psi.dot(tmp*psi).real();	
	return mean_value / norm_squared;
}

const double StrunzSimulatorReducedDeterministic::MAXIMUM_MAX_MIN_NORM_RATIO = 1.2;

// nested classes

StrunzSimulatorReducedDeterministic::WorkspaceHamiltonian::WorkspaceHamiltonian(const StrunzSimulatorReducedDeterministic& owner)
	: m_operators(owner.m_calc.data()), m_operators_next(owner.m_calc.data()), m_k1(owner.m_calc.data()), m_k2(owner.m_calc.data()), m_k3(owner.m_calc.data()), m_k4(owner.m_calc.data())
	, m_h1(owner.m_calc.nbrSites(), owner.m_calc.nbrSites()), m_h2(owner.m_calc.nbrSites(), owner.m_calc.nbrSites()), m_h3(owner.m_calc.nbrSites(), owner.m_calc.nbrSites()), m_h4(owner.m_calc.nbrSites(), owner.m_calc.nbrSites())
{
}

void StrunzSimulatorReducedDeterministic::WorkspaceHamiltonian::clear() 
{ 
	m_operators.clear(); 
}

StrunzSimulatorReducedDeterministic::WorkspaceState::WorkspaceState(const StrunzSimulatorReducedDeterministic& owner)
	: m_psi_k1(owner.m_calc.nbrSites()), m_psi_k2(owner.m_calc.nbrSites()), m_psi_k3(owner.m_calc.nbrSites()), m_psi_k4(owner.m_calc.nbrSites())
{
}

StrunzSimulatorReducedDeterministic::Workspace::Workspace(const StrunzSimulatorReducedDeterministic& owner)
	: m_hamiltonian_wksp(owner), m_state_wksp(owner)
{
}

