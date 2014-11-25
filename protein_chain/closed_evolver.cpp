#include "closed_evolver.h"
#include <Eigen/Eigenvalues>
#include <cassert>

ClosedEvolver::ClosedEvolver(const Eigen::MatrixXcd& hamiltonian)
  :  m_eigenvectors(hamiltonian.rows(), hamiltonian.cols()), m_eigenvalues(hamiltonian.rows())
{
  assert(hamiltonian.rows() == hamiltonian.cols());
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXcd> solver(hamiltonian);
  m_eigenvectors = solver.eigenvectors();
  m_eigenvalues = solver.eigenvalues();
}

void ClosedEvolver::evolve(const Eigen::VectorXcd& in, Eigen::VectorXcd& out, double t) const
{
  assert(in.size() == m_eigenvalues.size());
  assert(out.size() == m_eigenvalues.size());
  out.setZero();
  for (size_t evIdx = 0; evIdx < static_cast<unsigned int>(m_eigenvalues.size()); ++evIdx) {
    const std::complex<double> c(m_eigenvectors.col(evIdx).dot(in));
    const std::complex<double> phase = exp(std::complex<double>(0,-t) * m_eigenvalues[evIdx]);
    out += (c*phase)*m_eigenvectors.col(evIdx);
  }
}


