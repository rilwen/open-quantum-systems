#ifndef __EVOLVER_H
#define __EVOLVER_H

#include <Eigen/Core>
#include <vector>
#include "core.h"

class Evolver
{
public:
	PROTEIN_CHAIN_API virtual ~Evolver();
	virtual void evolve(const Eigen::VectorXcd& initState, const std::vector<double>& times, std::vector<Eigen::VectorXcd>& states) = 0;
	virtual void evolve(const Eigen::MatrixXcd& initStates, const std::vector<double>& times, std::vector<Eigen::MatrixXcd>& states) = 0;
};

#endif
