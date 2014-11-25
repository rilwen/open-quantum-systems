#include "evolver_stepping.h"
#include <cassert>
//#include <iostream>

EvolverStepping::EvolverStepping(double timeStep)
	: m_time_step(timeStep)
{
	assert( m_time_step > 0 );
}

void EvolverStepping::evolve(const Eigen::VectorXcd& initState, const std::vector<double>& times, std::vector<Eigen::VectorXcd>& states)
{
	evolveImpl<Eigen::VectorXcd>(initState, times, states);
}

void EvolverStepping::evolve(const Eigen::MatrixXcd& initStates, const std::vector<double>& times, std::vector<Eigen::MatrixXcd>& states)
{
	evolveImpl<Eigen::MatrixXcd>(initStates, times, states);
}

void EvolverStepping::evolve(const Eigen::VectorXcd& prev, double dt, Eigen::VectorXcd& next)
{
	evolveImpl<Eigen::VectorXcd>(prev, dt, next);
}

template <class State> void EvolverStepping::evolveImpl(const State& prev, const double delta_t, State& next)
{
	const unsigned int nbr_steps = static_cast<unsigned int>(ceil( delta_t / m_time_step ));
	//std::cout << "nbr_steps == " << nbr_steps << std::endl;
	next = prev;
	for (unsigned int k = 0; k < nbr_steps; ++k) {
		const double step_size = std::min( (k+1)*m_time_step, delta_t ) - k*m_time_step;
		assert( step_size >= 0 );
		assert( step_size <= m_time_step*(1+1E-6) );
		assert( k == nbr_steps - 1 || std::abs(step_size - m_time_step) < (1+1E-6)*m_time_step );
		step(step_size, next);
	}
}

template <class State> void EvolverStepping::evolveImpl(const State& initState, const std::vector<double>& times, std::vector<State>& states)
{
	const unsigned int len = times.size();
	if (len == 0) return;
	evolveImpl(initState, times[0], states[0]);
	for (unsigned int i = 1; i < len; ++i) {
		evolveImpl(states[i-1], times[i] - times[i-1], states[i]);
	}
}

void EvolverStepping::step(const double step_size, Eigen::MatrixXcd& states)
{
	Eigen::VectorXcd state(states.rows());
	for (unsigned int i = 0; i < static_cast<unsigned int>(states.cols()); ++i) {
		state = states.col(i);
		step(step_size, state);
		states.col(i) = state;
	}
}
