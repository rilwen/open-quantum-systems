#ifndef __STRUNZ_SIMULATOR_REDUCED_LINEAR_H
#define __STRUNZ_SIMULATOR_REDUCED_LINEAR_H

#include <boost/shared_ptr.hpp>
#include <vector>
#include <Eigen/Core>
#include "colored_gaussian_process_generator.h"
#include "strunz_simulator_reduced_deterministic.h"
#include "core.h"

class CorrelationFunctionDecomposable;

//! *Stochastic reduced-Strunz simulator (linear version without normalization)
class StrunzSimulatorReducedStochasticLinear
{
public:
	class Workspace
	{
	private:
		friend class StrunzSimulatorReducedStochasticLinear;
		Workspace(const StrunzSimulatorReducedStochasticLinear& owner);
		Eigen::MatrixXcd m_z_path;		
		StrunzSimulatorReducedDeterministic::WorkspaceState m_strunz_wksp_state;
		std::vector<Eigen::MatrixXcd> m_hamiltonians_times_minus_i_with_noise;
		boost::shared_ptr<const ColoredGaussianProcessGenerator> m_gauss_generator;
		boost::shared_ptr<ColoredGaussianProcessGenerator::Workspace> m_gauss_generator_wksp;		
		void reset();
	};
	QSD_API StrunzSimulatorReducedStochasticLinear(const std::vector<boost::shared_ptr<const CorrelationFunctionDecomposable> >& alphas, const Eigen::MatrixXcd& Hel, double dt, size_t nbr_steps);
	QSD_API StrunzSimulatorReducedStochasticLinear(boost::shared_ptr<const CorrelationFunctionDecomposable> alpha, size_t nbr_sites, const Eigen::MatrixXcd& Hel, double dt, size_t nbr_steps);
	QSD_API Workspace workspace() const;
	QSD_API void simulate(Workspace& wksp, const Eigen::VectorXcd& initial, std::vector<Eigen::VectorXcd>& trajectory) const;
private:
	void init();
	void simulate_z(Workspace& wksp) const;
	void calculate_hamiltonians_with_noise(Workspace& wksp) const;
private:	
	StrunzSimulatorReducedDeterministic m_strunz;
	double m_dt;
	size_t m_nbr_steps;
	std::vector<Eigen::MatrixXcd> m_hamiltonians_times_minus_i_no_noise;
	bool m_common_alpha;
};

#endif // __STRUNZ_SIMULATOR_REDUCED_LINEAR_H
