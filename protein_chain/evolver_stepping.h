#ifndef __EVOLVER_STEPPING_H
#define __EVOLVER_STEPPING_H

#include "evolver.h"
#include "core.h"
#include <Eigen/Core>

class EvolverStepping: public Evolver
{
public:
	PROTEIN_CHAIN_API EvolverStepping(double timeStep);
	PROTEIN_CHAIN_API void evolve(const Eigen::VectorXcd& initState, const std::vector<double>& times, std::vector<Eigen::VectorXcd>& states);
	PROTEIN_CHAIN_API void evolve(const Eigen::MatrixXcd& initStates, const std::vector<double>& times, std::vector<Eigen::MatrixXcd>& states);
	PROTEIN_CHAIN_API void evolve(const Eigen::VectorXcd& prev, double dt, Eigen::VectorXcd& next);	
	template <class Acc> void evolve(const Eigen::VectorXcd& prev, double dt, Eigen::VectorXcd& next, Acc& accumulator);	
	PROTEIN_CHAIN_API virtual void step(double step_size, Eigen::VectorXcd& state) = 0; // single step of the method
	double timeStep() const { return m_time_step; }
private:
	template <class State> void evolveImpl(const State& initState, const std::vector<double>& times, std::vector<State>& states);
	template <class State> void evolveImpl(const State& prev, double delta_t,  State& next);
	template <class State, class Acc> void evolveImpl(const State& prev, double delta_t,  State& next, Acc& accumulator);
	PROTEIN_CHAIN_API virtual void step(double step_size, Eigen::MatrixXcd& states); // single step of the method for a matrix of vectors (default implementation provided)
private:
	double m_time_step;
};

template <class State, class Acc> void EvolverStepping::evolveImpl(const State& prev, const double delta_t, State& next, Acc& accumulator)
{
	const unsigned int nbr_steps = static_cast<unsigned int>(ceil( delta_t / m_time_step ));
	//std::cout << "nbr_steps == " << nbr_steps << std::endl;
	accumulator(0, prev);
	next = prev;
	for (unsigned int k = 0; k < nbr_steps; ++k) {
		const double new_t = std::min( (k+1)*m_time_step, delta_t );
		const double step_size = new_t - k*m_time_step;
		assert( step_size >= 0 );
		assert( step_size <= m_time_step*(1+1E-6) );
		assert( k == nbr_steps - 1 || std::abs(step_size - m_time_step) < (1+1E-6)*m_time_step );
		step(step_size, next);
		accumulator(new_t, next);
	}
}

template <class Acc> void EvolverStepping::evolve(const Eigen::VectorXcd& prev, double dt, Eigen::VectorXcd& next, Acc& accumulator)
{
	evolveImpl(prev, dt, next, accumulator);
}

#endif // __EVOLVER_STEPPING_H
