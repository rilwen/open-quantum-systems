#ifndef __STRUNZ_SIMULATOR_REDUCED_NONLINEAR_H
#define __STRUNZ_SIMULATOR_REDUCED_NONLINEAR_H

#include <boost/shared_ptr.hpp>
#include <vector>
#include <Eigen/Core>
#include <math/Jagged2DArray.h>
#include "colored_gaussian_process_generator.h"
#include "strunz_simulator_reduced_deterministic.h"
#include "core.h"

class CorrelationFunctionDecomposable;

//! *Stochastic reduced-Strunz simulator (nonlinear version with normalization)
class StrunzSimulatorReducedStochasticNonlinear
{
public:
	//! Psi state + elementary drift integrals
	class State
	{
	public:
		QSD_API State() {}
		QSD_API State& operator+=(const State& rhs);
		QSD_API State& operator-=(const State& rhs);
		QSD_API State& operator*=(const std::complex<double>& x);
		QSD_API State& operator/=(const std::complex<double>& x);
		QSD_API double norm() const { return m_psi.norm(); }
	private:
		friend class StrunzSimulatorReducedStochasticNonlinear;		
		State(const StrunzSimulatorReducedStochasticNonlinear& owner);
		Eigen::VectorXcd m_psi;
		rql::math::Jagged2DArray<std::complex<double> > m_drift_integrals;
		

		//! set new psi, zero the drift integrals and initialize the interaction energy at time 0
		void reset(const Eigen::VectorXcd& initial, const StrunzSimulatorReducedStochasticNonlinear& owner);
	};

	class Workspace
	{
	private:
		friend class StrunzSimulatorReducedStochasticNonlinear;
		Workspace(const StrunzSimulatorReducedStochasticNonlinear& owner);
		const StrunzSimulatorReducedStochasticNonlinear& m_owner;
		Eigen::MatrixXcd m_z_path;		
		boost::shared_ptr<const ColoredGaussianProcessGenerator> m_gauss_generator;
		boost::shared_ptr<ColoredGaussianProcessGenerator::Workspace> m_gauss_generator_wksp;
		Eigen::MatrixXcd m_nonlinear_hamiltonian_times_i;
		Eigen::MatrixXcd m_nonlinear_hamiltonian;
		Eigen::MatrixXcd m_interaction_hamiltonian;
		Eigen::VectorXcd m_tmp_vec;
		State m_state;
		void reset();
	};

	class Functor
	{
	public:				
		void operator()(double t, const State& y, State& dy);
	private:
		Functor(Workspace& wksp);
		friend class StrunzSimulatorReducedStochasticNonlinear;

		size_t m_nbr_sites;
		// not safe!
		Workspace& m_wksp;
		const CorrelationFunctionsDecomposition& m_alpha_decomp;
	};

	QSD_API StrunzSimulatorReducedStochasticNonlinear(const std::vector<boost::shared_ptr<const CorrelationFunctionDecomposable> >& alphas, const Eigen::MatrixXcd& Hel, double dt, size_t nbr_steps);
	QSD_API StrunzSimulatorReducedStochasticNonlinear(boost::shared_ptr<const CorrelationFunctionDecomposable> alpha, size_t nbr_sites, const Eigen::MatrixXcd& Hel, double dt, size_t nbr_steps);
	QSD_API Workspace workspace() const;
	QSD_API void simulate(Workspace& wksp, const Eigen::VectorXcd& initial, std::vector<Eigen::VectorXcd>& trajectory, std::vector<double>* interaction_energies = 0) const;
private:
	void init(double dt, size_t nbr_steps);
	void simulate_z(Workspace& wksp) const;
	void calc_time_interpolation_weights_and_indices(double t, double dt, double& prev_operator_weight, size_t& prev_operator_idx, double& next_operator_weight, size_t& next_operator_idx) const;
	//! contains intermediate results calculated during interpolation
	struct InterpolationParams
	{
		size_t noise_idx;
		size_t prev_ham_idx;
		size_t next_ham_idx;
		double weight_prev;
		double weight_next; 
		size_t prev_operator_start_idx;
		size_t next_operator_start_idx;
	};
	InterpolationParams interpolate_calculator_outputs(Workspace& wksp, double t, const std::vector<Eigen::MatrixXcd>& outputs, Eigen::MatrixXcd& result) const;
	double calc_interaction_energy(Workspace& wksp, double t, const State& y) const;
private:	
	StrunzSimulatorReducedDeterministic m_strunz;
	double m_dt;
	size_t m_nbr_steps;
	size_t m_nbr_hamiltonians;
	std::vector<Eigen::MatrixXcd> m_hamiltonians_times_minus_i_no_noise;
	std::vector<Eigen::MatrixXcd> m_interaction_hamiltonians_no_noise;
	std::vector<Eigen::MatrixXcd> m_operators;
	bool m_common_alpha;
};

#endif // __STRUNZ_SIMULATOR_REDUCED_NONLINEAR_H
