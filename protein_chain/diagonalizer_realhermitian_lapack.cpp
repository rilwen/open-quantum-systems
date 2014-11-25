#include "diagonalizer_realhermitian_lapack.h"
#include <sstream>
#include <stdexcept>

extern "C" {
	void dsyevd_(const char* jobz, const char* uplo, int* n, double* a, int* lda, double* w, double* work, int* lwork, int* iwork, int* liwork, int* info);
}

DiagonalizerRealHermitianLapack::DiagonalizerRealHermitianLapack(unsigned int dim)
	: m_dim(dim), m_lwork(-1), m_liwork(1), m_work(new double[1]), m_iwork(new int[1])
{
	m_info = 0;
	boost::scoped_array<double> w(new double[dim]);
	dsyevd_("V", "U", &m_dim, 0, &m_dim, w.get(), m_work.get(), &m_lwork, m_iwork.get(), &m_liwork, &m_info);
	if (m_info) {
		std::stringstream ss;
		ss << "DSYEVD call failed in constructor with info == " << m_info;
		throw std::runtime_error(ss.str());
	}
	m_lwork = static_cast<int>(m_work[0]);
	m_liwork = m_iwork[0];
	m_work.reset(new double[m_lwork]);
	m_iwork.reset(new int[m_liwork]);
}

void DiagonalizerRealHermitianLapack::diagonalize(const Eigen::MatrixXd& m, Eigen::MatrixXd& eigenvectors, Eigen::VectorXd& eigenvalues)
{
	eigenvectors = m;
	dsyevd_("V", "U", &m_dim, eigenvectors.data(), &m_dim, eigenvalues.data(), m_work.get(), &m_lwork, m_iwork.get(), &m_liwork, &m_info);
	if (m_info) {
		std::stringstream ss;
		ss << "DSYEVD call failed in constructor with diagonalize == " << m_info;
		throw std::runtime_error(ss.str());
	}
}
