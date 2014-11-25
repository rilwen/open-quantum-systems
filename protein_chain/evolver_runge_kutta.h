#ifndef __EVOLVER_RUNGE_KUTTA_H
#define __EVOLVER_RUNGE_KUTTA_H

#include "evolver_stepping.h"
#include "core.h"
#include <Eigen/Core>
#include <boost/array.hpp>

class EvolverRungeKutta: public EvolverStepping
{
public:
	PROTEIN_CHAIN_API EvolverRungeKutta(const Eigen::MatrixXcd& hamiltonian, double timeStep);
	PROTEIN_CHAIN_API_DEBUG void step(double step_size, Eigen::VectorXcd& state); // single step of the RK method
private:	
	PROTEIN_CHAIN_API_DEBUG void step(double step_size, Eigen::MatrixXcd& states); // single step of the RK method
private:
	static const unsigned int m_order = 4;
	Eigen::MatrixXcd m_M; // M = -i*H
	boost::array<Eigen::VectorXcd,m_order> m_k;
	static const boost::array<double,4> m_b;
	static const boost::array<double,4> m_c;
	static const boost::array<double,1> m_a1;
	static const boost::array<double,2> m_a2;
	static const boost::array<double,3> m_a3;
};

#endif // __EVOLVER_RUNGE_KUTTA_H

