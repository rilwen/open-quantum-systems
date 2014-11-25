#ifndef __EVOLVERNORMAL_H
#define __EVOLVERNORMAL_H

#include <Eigen/Core>
#include <vector>
#include "evolver.h"
#include <Eigen/Eigenvalues>

class EvolverNormal: public Evolver
{
public:
	EvolverNormal(const Eigen::ComplexEigenSolver<MatrixXcd>& diagonalizer);
	virtual void evolve(const Eigen::VectorXcd& initState, const std::vector<double>& times, std::vector<Eigen::VectorXcd>& states);
	virtual void evolve(const Eigen::MatrixXcd& initStates, const std::vector<double>& times, std::vector<Eigen::MatrixXcd>& states);
private:
	Eigen::ComplexEigenSolver<MatrixXcd> m_diagonalizer;
};

#endif
