#ifndef __STRUNZ_SIMULATOR_GRID_H
#define __STRUNZ_SIMULATOR_GRID_H

#include <boost/shared_ptr.hpp>
#include <vector>
#include <Eigen/Core>
#include "strunz_calculator_grid_explicit.h"
#include "colored_gaussian_process_generator.h"
#include "core.h"

class CorrelationFunction;


class StrunzSimulatorGrid
{
public:
	class Workspace
	{
	private:
		friend class StrunzSimulatorGrid;
		Workspace(const StrunzSimulatorGrid& owner);		
		Eigen::MatrixXcd m_Heff_with_noise;
		Eigen::MatrixXcd m_z_path;
		boost::shared_ptr<const ColoredGaussianProcessGenerator> m_gauss_generator;
		boost::shared_ptr<ColoredGaussianProcessGenerator::Workspace> m_gauss_generator_wksp;
	};
	//! @param[in] stochastic Perform stochastic simulation
	QSD_API_DEBUG StrunzSimulatorGrid(const std::vector<boost::shared_ptr<const CorrelationFunction> >& alpha, const Eigen::MatrixXcd& Hel, double dt, size_t nbrTimeSteps, size_t nbrCalcStepsPerStep, bool stochastic);
	//! @param[in] stochastic Perform stochastic simulation
	QSD_API_DEBUG StrunzSimulatorGrid(boost::shared_ptr<const CorrelationFunction> alpha, size_t dim, const Eigen::MatrixXcd& Hel, double dt, size_t nbrTimeSteps, size_t nbrCalcStepsPerStep, bool stochastic);
	QSD_API_DEBUG Workspace workspace() const;
	size_t dim() const { return m_nbr_sites; }
	QSD_API_DEBUG void simulateTrajectory(Workspace& wksp, const Eigen::VectorXcd& initial, std::vector<Eigen::VectorXcd>& trajectory) const;
	size_t nbrSteps() const { return m_nbr_steps; }
	double dt() const { return m_dt; }
	const Eigen::MatrixXcd& effective_hamiltonian_without_noise(size_t i) const { return m_Heff[i]; }
private:
	void calculate_effective_hamiltonians_without_noise();
	//! @param[in] dt Time step for the SDE discretization
	void calculate_effective_hamiltonian_with_noise(Workspace& wksp, size_t sim_idx, double dt) const;
	double taylor_evolver_time_step(const Eigen::MatrixXcd& effective_hamiltonian, double base_dt) const;
	//! @param[in] dt Time step for the SDE discretization
	const Eigen::MatrixXcd& hamiltonian_for_evolution(Workspace& wksp, size_t sim_idx, double dt) const;
	void simulate_z(Workspace& wksp) const;
	void reset_z(Workspace& wksp) const { wksp.m_z_path.setZero(); }
	StrunzSimulatorGrid& operator=(const StrunzSimulatorGrid&); // not implemented
private:
	std::vector<boost::shared_ptr<const CorrelationFunction> > m_alpha;
	const Eigen::MatrixXcd m_Hel;
	const size_t m_nbr_sites;
	const double m_dt;	
	const size_t m_nbr_steps;
	const size_t m_nbr_calc_steps_per_step;
	const double m_calc_dt;
	const size_t m_nbr_calc_steps;
	const bool m_stochastic;	
	const bool m_common_alpha;
	const double m_Hel_norm;
	std::vector<Eigen::MatrixXcd> m_Heff;
};

#endif // __STRUNZ_SIMULATOR_GRID_H
