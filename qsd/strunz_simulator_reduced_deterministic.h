#ifndef __QSD_STRUNZ_ABSORPTION_SPECTRUM_REDUCED_H
#define __QSD_STRUNZ_ABSORPTION_SPECTRUM_REDUCED_H

#include <utility>
#include <vector>
#include <boost/shared_ptr.hpp>
#include <boost/array.hpp>
#include <Eigen/Core>
#include "strunz_calculator_reduced.h"
#include "core.h"

class CorrelationFunctionDecomposable;

//! Simulates psi and effective hamiltonians without Gaussian noise, integrating the evolution equations using 4th order Runge-Kutta method
class StrunzSimulatorReducedDeterministic
{
private:
	static const unsigned int MAX_ALLOWED_DT_DIVISIONS = 10;
public:
	//! Workspace used to evolve the effective hamiltonian
	class WorkspaceHamiltonian
	{
	private:
		WorkspaceHamiltonian(const StrunzSimulatorReducedDeterministic& owner);
		friend class StrunzSimulatorReducedDeterministic;
		StrunzCalculatorReduced::Data m_operators;
		StrunzCalculatorReduced::Data m_operators_next;
		StrunzCalculatorReduced::Data m_k1;	
		StrunzCalculatorReduced::Data m_k2;
		StrunzCalculatorReduced::Data m_k3;
		StrunzCalculatorReduced::Data m_k4;
		Eigen::MatrixXcd m_h1;
		Eigen::MatrixXcd m_h2;
		Eigen::MatrixXcd m_h3;
		Eigen::MatrixXcd m_h4;
		void clear();
	};

	//! Workspace used to evolve the state
	class WorkspaceState
	{
	private:
		WorkspaceState(const StrunzSimulatorReducedDeterministic& owner);
		friend class StrunzSimulatorReducedDeterministic;
		Eigen::VectorXcd m_psi_k1;
		Eigen::VectorXcd m_psi_k2;
		Eigen::VectorXcd m_psi_k3;
		Eigen::VectorXcd m_psi_k4;
		Eigen::MatrixXcd m_tmp;
		boost::array<Eigen::MatrixXcd, MAX_ALLOWED_DT_DIVISIONS> m_mid_h;
		boost::array<Eigen::VectorXcd, MAX_ALLOWED_DT_DIVISIONS> m_mid_psi;
	};

	//! Workspace combining both hamiltonian and state evolution
	class Workspace
	{
	private:
		Workspace(const StrunzSimulatorReducedDeterministic& owner);
		friend class StrunzSimulatorReducedDeterministic;
		WorkspaceHamiltonian& hamiltonian_workspace() { return m_hamiltonian_wksp; }		
		const WorkspaceHamiltonian& hamiltonian_workspace() const { return m_hamiltonian_wksp; }
		WorkspaceState& state_workspace() { return m_state_wksp; }
		const WorkspaceState& state_workspace() const { return m_state_wksp; }
		void clear() { m_hamiltonian_wksp.clear(); }
		WorkspaceHamiltonian m_hamiltonian_wksp;
		WorkspaceState m_state_wksp;
	};

	//! Acceptor for diagnostic data about the simulation
	class Diagnostics
	{
	public:
		typedef std::pair<double,double> data_type;
		typedef std::vector<data_type> data_store_type;
		void add_psi_norm(double t, double norm) { m_psi_norms.push_back(data_type(t, norm)); }
		void add_trace_deviation(double t, double deviation) { m_trace_deviations.push_back(data_type(t, deviation)); }
		void add_hamiltonian_norm(double t, double norm) { m_hamiltonian_norms.push_back(data_type(t, norm)); }
		const data_store_type& psi_norms() const { return m_psi_norms; }
		const data_store_type& trace_deviations() const { return m_trace_deviations; }
		const data_store_type& hamiltonian_norms() const { return m_hamiltonian_norms; }
	private:
		data_store_type m_psi_norms;
		data_store_type m_trace_deviations;
		data_store_type m_hamiltonian_norms;
	};

	QSD_API StrunzSimulatorReducedDeterministic(const std::vector<boost::shared_ptr<const CorrelationFunctionDecomposable> >& alphas, const Eigen::MatrixXcd& Hel, bool correct_psi_norm);
	QSD_API Workspace workspace() const { return Workspace(*this); }
	QSD_API WorkspaceHamiltonian workspace_hamiltonian() const { return WorkspaceHamiltonian(*this); }
	QSD_API WorkspaceState workspace_state() const { return WorkspaceState(*this); }
	QSD_API void simulate(Workspace& wksp, const Eigen::VectorXcd& initial, std::vector<Eigen::VectorXcd>& trajectory, double dt, size_t nbr_steps, Diagnostics* diagnostics = 0) const;
	//! fill scalar_products with <psi(0)|psi(t)>
	QSD_API void simulate(Workspace& wksp, const Eigen::VectorXcd& initial, std::vector<std::complex<double> >& scalar_products, double dt, size_t nbr_steps) const;
	//! Simulate effective hamiltonian (multiplied by -I) non-adaptively (single Runge-Kutta step of length dt).
	//! @param[in] all If true, save all 4 hamiltonians per step and the last one (4*nbr_steps + 1 in total). If false, save only the first one for each step and, the last one (nbr_step + 1 in total)
	QSD_API void simulate_hamiltonians_times_minus_i_nonadaptive(WorkspaceHamiltonian& wksp, std::vector<Eigen::MatrixXcd>& hamiltonians_times_minus_i, double dt, size_t nbr_steps, bool all) const;
	//! Simulate effective hamiltonians (multiplied by -I) adaptive
	QSD_API void simulate_hamiltonians_times_minus_i_adaptive(WorkspaceHamiltonian& wksp, std::vector<Eigen::MatrixXcd>& hamiltonians_times_minus_i, double dt, size_t nbr_steps, std::vector<Eigen::MatrixXcd>* interaction_hamiltonians = 0) const;
	//! Perform the same simulation as simulate_hamiltonians_times_minus_i_adaptive, but at the end of each step call listener(const_cast<const WorkspaceHamiltonian& >(wksp), idx) where idx = 0, ..., nbr_steps - 1 is the step index.
	//! Listener should only read the workspace data
	template <class Listener> void simulate_hamiltonians_times_minus_i_adaptive(WorkspaceHamiltonian& wksp, std::vector<Eigen::MatrixXcd>& hamiltonians_times_minus_i, double dt, size_t nbr_steps, Listener& listener, std::vector<Eigen::MatrixXcd>* interaction_hamiltonians = 0) const;
	//! Simulate state using pre-calculated hamiltonians (needs all 4 hamiltonians for 1 RK step), multiplied by -I
	//! @param[in] hamiltonians_times_minus_i Vector of hamiltonians, as returned by simulate_hamiltonians_times_minus_i method with all=TRUE (last hamiltonian can be truncated); hamiltonians.size() >= 4*nbr_steps
	QSD_API void simulate_psi(WorkspaceState& wksp, const std::vector<Eigen::MatrixXcd>& hamiltonians_times_minus_i, const Eigen::VectorXcd& initial, std::vector<Eigen::VectorXcd>& trajectory, double dt, size_t nbr_steps) const;
	QSD_API size_t nbr_sites() const { return m_calc.nbrSites(); }
	QSD_API const std::vector<boost::shared_ptr<const CorrelationFunctionDecomposable> >& alphas() const { return m_alphas; }
	QSD_API void calculate_operator(const WorkspaceHamiltonian& wksp, size_t site_idx, Eigen::MatrixXcd& m) const { m_calc.calculate_operator(wksp.m_operators, site_idx, m); }
	QSD_API const CorrelationFunctionsDecomposition& correlation_functions_decomposition() const { return m_calc.correlation_functions_decomposition(); }
private:
	//! Perform 1 step of 4th order Runge-Kutta effective hamiltonian evolution
	//! @param force Always update the state
	//! @return TRUE if the convergence criteria were satisfied (and wksp.m_operators contains the new set of states), FALSE if the convergence criteria were not satisfied (and wksp.m_operators was not updated, unless force = TRUE)	
	bool step_hamiltonian_impl(WorkspaceHamiltonian& wksp, double t, double dt, bool force) const;
	//! Perform 1 step of 4th order Runge-Kutta effective hamiltonian evolution (recursive function which changes dt adaptively to ensure better accuracy)
	void step_hamiltonian_adaptive(WorkspaceHamiltonian& wksp, double t, double dt, unsigned int nbr_allowed_dt_divisions) const;
	//! Perform 1 step of 4th order Runge-Kutta state evolution, assuming that wksp_ham contains the data from the corresponding step of the hamiltonian evolution
	void step_psi_nonadaptive(const WorkspaceHamiltonian& wksp_ham, WorkspaceState& wksp_st, const Eigen::VectorXcd& prev, Eigen::VectorXcd& next, double dt) const
	{
		step_psi_impl(wksp_ham.m_h1, wksp_ham.m_h2, wksp_ham.m_h3, wksp_ham.m_h4, wksp_st, prev, next, dt, true); // important: don't use anything else from wksp_ham then m_h1,m_h2,m_h3,m_h4
	}
	//! Perform 1 step of 4th order Runge-Kutta state evolution ADAPTIVELY, assuming that wksp_ham contains the data from the corresponding step of the hamiltonian evolution
	void step_psi_adaptive(const WorkspaceHamiltonian& wksp_ham, WorkspaceState& wksp_st, const Eigen::VectorXcd& prev, Eigen::VectorXcd& next, double dt) const
	{
		step_psi_adaptive(wksp_ham.m_h1, wksp_ham.m_h2, wksp_ham.m_h3, wksp_ham.m_h4, wksp_st, prev, next, dt, MAX_ALLOWED_DT_DIVISIONS);
	}
	//! Perform 1 step of 4th order Runge-Kutta state evolution, using the precalculated effective hamiltonians (multiplied by -I) from the corresponding step of the hamiltonian evolution
	//! Corrects the norm of the evolved state
	//! @return TRUE if the convergence criteria were satisfied (and next contains the new psi), FALSE if the convergence criteria were not satisfied (and next was not updated, unless force = TRUE)
	bool step_psi_impl(const Eigen::MatrixXcd& h1, const Eigen::MatrixXcd& h2, const Eigen::MatrixXcd& h3, const Eigen::MatrixXcd& h4, WorkspaceState& wksp_st, const Eigen::VectorXcd& prev, Eigen::VectorXcd& next, double dt, bool force) const;
	//! Perform 1 step (adaptively decreasing dt if necessary) of 4th order Runge-Kutta state evolution, using the precalculated effective hamiltonians (multiplied by -I) from the corresponding step of the hamiltonian evolution
	//! Corrects the norm of the evolved state
	void step_psi_adaptive(const Eigen::MatrixXcd& h1, const Eigen::MatrixXcd& h2, const Eigen::MatrixXcd& h3, const Eigen::MatrixXcd& h4, WorkspaceState& wksp_st, const Eigen::VectorXcd& prev, Eigen::VectorXcd& next, double dt, unsigned int nbr_allowed_dt_divisions) const;
	//! Evolve both hamiltonian and state
	void step_both(Workspace& wksp, const Eigen::VectorXcd& prev, Eigen::VectorXcd& next, double t, double dt) const
	{
		step_both(wksp.hamiltonian_workspace(), wksp.state_workspace(), prev, next, t, dt);
	}
	//! Evolve both hamiltonian and state
	void step_both(WorkspaceHamiltonian& wksp_ham, WorkspaceState& wksp_st, const Eigen::VectorXcd& prev, Eigen::VectorXcd& next, double t, double dt) const
	{
		step_both_adaptive(wksp_ham, wksp_st, prev, next, t, dt, MAX_ALLOWED_DT_DIVISIONS);
	}
	//! Evolve both hamiltonian and state (recursive function)
	//! @param nbr_allowed_dt_divisions Number of times it is still allowed to cut the time step in half
	void step_both_adaptive(WorkspaceHamiltonian& wksp_ham, WorkspaceState& wksp_st, const Eigen::VectorXcd& prev, Eigen::VectorXcd& next, double t, double dt, unsigned int nbr_allowed_dt_divisions) const;
	//! Calculate time derivative of ln<Psi|Psi>. Does not depend on norm of psi
	double calculate_log_norm_derivative(const Eigen::VectorXcd& psi, const Eigen::MatrixXcd& h_eff_times_minus_i, Eigen::MatrixXcd& tmp) const;
	//! Calculate time derivative of ln<Psi|Psi>. Does not depend on norm of psi
	double calculate_log_norm_derivative(WorkspaceState& wksp, const Eigen::VectorXcd& psi, const Eigen::MatrixXcd& h_eff_times_minus_i) const
	{
		return calculate_log_norm_derivative(psi, h_eff_times_minus_i, wksp.m_tmp);
	}
	static bool norms_diverged(double min_norm, double max_norm);
private:	
	StrunzCalculatorReduced m_calc;		
	std::vector<boost::shared_ptr<const CorrelationFunctionDecomposable> > m_alphas;
	bool m_correct_psi_norm;
	static const double MAXIMUM_MAX_MIN_NORM_RATIO;
};

template <class Listener> void StrunzSimulatorReducedDeterministic::simulate_hamiltonians_times_minus_i_adaptive(WorkspaceHamiltonian& wksp, std::vector<Eigen::MatrixXcd>& hamiltonians_times_minus_i, const double dt, const size_t nbr_steps, Listener& listener, std::vector<Eigen::MatrixXcd>* interaction_hamiltonians) const
{
	wksp.clear();
	hamiltonians_times_minus_i.resize(nbr_steps + 1);
	for (size_t i = 0; i < nbr_steps; ++i) {
		m_calc.calculate_effective_hamiltonian_times_minus_i(wksp.m_operators, hamiltonians_times_minus_i[i]);
		if (interaction_hamiltonians) {
			m_calc.calculate_interaction_hamiltonian(wksp.m_operators, (*interaction_hamiltonians)[i]);
		}
		step_hamiltonian_adaptive(wksp, i*dt, dt, MAX_ALLOWED_DT_DIVISIONS);
		listener(const_cast<const WorkspaceHamiltonian& >(wksp), i);
	}
	m_calc.calculate_effective_hamiltonian_times_minus_i(wksp.m_operators, hamiltonians_times_minus_i.back());
	if (interaction_hamiltonians) {
		m_calc.calculate_interaction_hamiltonian(wksp.m_operators, interaction_hamiltonians->back());
	}
}


#endif // __QSD_STRUNZ_ABSORPTION_SPECTRUM_REDUCED_H
