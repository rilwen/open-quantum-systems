#ifndef __CLOSED_EVOLVER_H
#define __CLOSED_EVOLVER_H

#include <Eigen/Core>

class ClosedEvolver
{
public:
  ClosedEvolver(const Eigen::MatrixXcd& hamiltonian);
  //! @param[in] t evolution time
  void evolve(const Eigen::VectorXcd& in, Eigen::VectorXcd& out, double t) const;
  size_t dimension() const { return m_eigenvalues.size(); }
private:
  Eigen::MatrixXcd m_eigenvectors;
  Eigen::VectorXd m_eigenvalues;
};

#endif // __CLOSED_EVOLVER_H
