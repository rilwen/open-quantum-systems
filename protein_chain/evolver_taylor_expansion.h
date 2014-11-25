#ifndef __EVOLVER_TAYLOR_EXPANSION_H
#define __EVOLVER_TAYLOR_EXPANSION_H

#include "evolver_stepping.h"
#include <Eigen/Core>
#include <boost/array.hpp>
#include <cassert>
#include <cmath>

template <unsigned int Order = 4u>
class EvolverTaylorExpansion: public EvolverStepping
{
public:
	EvolverTaylorExpansion(const Eigen::MatrixXcd& hamiltonian, double timeStep);
	void step(double step_size, Eigen::VectorXcd& state);
private:
	Eigen::MatrixXcd m_standard_step;
	boost::array<Eigen::MatrixXcd,Order+1> m_iHpow;
	Eigen::VectorXcd m_work;
};

template <unsigned int Order> EvolverTaylorExpansion<Order>::EvolverTaylorExpansion(const Eigen::MatrixXcd& hamiltonian, double timeStep)
	: EvolverStepping(timeStep), m_work(hamiltonian.rows())
{
	assert( hamiltonian.rows() == hamiltonian.cols() );
	const Eigen::MatrixXcd minusIH(hamiltonian*std::complex<double>(0.0,-1.0));
	m_standard_step = Eigen::MatrixXcd::Identity(hamiltonian.rows(), hamiltonian.rows());
	m_iHpow[0] = m_standard_step;
	double tpowa = timeStep;
	for (unsigned int k = 1; k <= Order; ++k) {
		m_iHpow[k] = m_iHpow[k-1] * minusIH;
		m_iHpow[k] /= k;
		m_standard_step += m_iHpow[k] * tpowa;
		tpowa *= timeStep;
	}
}

template <unsigned int Order> void EvolverTaylorExpansion<Order>::step(double step_size, Eigen::VectorXcd& state)
{	
	if (step_size == timeStep()) {
		state = (m_standard_step*state);
	} else {
		m_work = state;
		double hpow = 1;
		for (unsigned int k = 1; k <= Order; ++k) {
			hpow *= step_size;
			state.noalias() += (m_iHpow[k]*m_work)*hpow;
		}
	}
}

#endif // __EVOLVER_TAYLOR_EXPANSION_H
