#ifndef __RHO_ANALYZER_HPP
#define __RHO_ANALYZER_HPP

#include <Eigen/Core>
#include "core.h"

class RhoAnalyzer
{
public:
	PROTEIN_CHAIN_API RhoAnalyzer(const Eigen::MatrixXcd& rho, bool wrap = false, double epsilon = 1E-10);
	PROTEIN_CHAIN_API double vonNeumannEntropy() const;
	double vonNeumannEntropy(double epsilon) const;
	PROTEIN_CHAIN_API double coherence() const;
	const Eigen::VectorXd& eigenvalues() const { return m_eigenvalues; }
	const Eigen::VectorXd& original_eigenvalues() const { return m_eigenvalues_original; }
	const Eigen::MatrixXcd& eigenvectors() const { return m_eigenvectors; }
	void correct_rho(Eigen::MatrixXcd& rho) const;
	void original_rho(Eigen::MatrixXcd& rho) const;
	bool were_bad_evals_found() const { return m_bad_evals_found; }
private:
	bool m_wrap;
	Eigen::MatrixXcd m_eigenvectors;
	Eigen::VectorXd m_eigenvalues_original;
	Eigen::VectorXd m_eigenvalues;
	size_t m_dim;
	double m_epsilon;
	bool m_bad_evals_found;
};

#endif // __RHO_ANALYZER_HPP
