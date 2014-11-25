#include "rho_analyzer.h"
#include <Eigen/Eigenvalues>
#include <cmath>
#include <cassert>
#include <stdexcept>
#include <sstream>
#include <iostream>
#include <boost/math/special_functions/fpclassify.hpp>

static double fix_eigenvalue(double l)
{	
	if (l < 0) {
		return fix_eigenvalue(-l);
	} else if (l > 1) {
		return fix_eigenvalue(2 - l);
	} else {
		return l;
	}
}

static inline double is_eval_ok(double l, double epsilon) {
	return (l >= -epsilon && l <= 1 + epsilon);
}

RhoAnalyzer::RhoAnalyzer(const Eigen::MatrixXcd& rho, const bool wrap, const double epsilon)
	: m_wrap(wrap), m_dim(rho.rows()), m_epsilon(epsilon)
{
	assert(rho.rows() == rho.cols());
	Eigen::SelfAdjointEigenSolver<Eigen::MatrixXcd> solver(rho);
	m_eigenvectors = solver.eigenvectors();
	m_eigenvalues_original = solver.eigenvalues();
	if (wrap) {
		double trace = 0;
		m_eigenvalues.resize(m_eigenvalues_original.size());
		for (size_t i = 0; i < static_cast<size_t>(rho.rows()); ++i) {
			m_eigenvalues[i] = fix_eigenvalue(m_eigenvalues_original[i]);
			trace += m_eigenvalues[i];
		}
		m_eigenvalues /= trace;
	} else {
		m_eigenvalues = m_eigenvalues_original;
	}
	m_bad_evals_found = false;
	for (size_t i = 0; i < static_cast<size_t>(rho.rows()); ++i) {
		m_bad_evals_found |= (! is_eval_ok(m_eigenvalues_original[i], m_epsilon));
	}
}

double RhoAnalyzer::vonNeumannEntropy() const
{
	return vonNeumannEntropy(m_epsilon);
}

double RhoAnalyzer::vonNeumannEntropy(double epsilon) const
{
	double entropy = 0;
	const double tolerance = std::max(1.0, m_eigenvalues.array().abs().maxCoeff()) * epsilon;	
	for (size_t i = 0; i < m_dim; ++i) {
		const double p = m_eigenvalues[i];
		if (boost::math::isnan(p))
			throw std::runtime_error("Eigenvalue is NaN");
		if (p > 0) {
			entropy -= p * log(p);
		} else if (p < -tolerance) {
			std::stringstream ss;
			ss << "Negative rho eigenvalue: " << p << " exceeding tolerance " << tolerance;
			throw std::runtime_error(ss.str().c_str());
		}
	}
	return entropy;
}

double RhoAnalyzer::coherence() const
{
	double coh = 0;
	for (size_t i = 0; i < m_dim; ++i) {
		const double p = m_eigenvalues[i];
		coh += p*p;
	}
	return coh;
}

void RhoAnalyzer::correct_rho(Eigen::MatrixXcd& rho) const
{
	rho.setZero(m_dim, m_dim);
	for (size_t i = 0; i < m_dim; ++i) {
		const double corr_l = m_wrap ? m_eigenvalues[i] : fix_eigenvalue(m_eigenvalues[i]);
		rho += corr_l * m_eigenvectors.col(i) * m_eigenvectors.col(i).adjoint();
	}
	rho /= rho.trace();
}

void RhoAnalyzer::original_rho(Eigen::MatrixXcd& rho) const
{
	rho.setZero(m_dim, m_dim);
	for (size_t i = 0; i < m_dim; ++i) {
		rho += m_eigenvalues_original[i] * m_eigenvectors.col(i) * m_eigenvectors.col(i).adjoint();
	}
	rho /= rho.trace();
}

