#ifndef __DIAGONALIZER_REALHERMITIAN_LAPACK_H
#define __DIAGONALIZER_REALHERMITIAN_LAPACK_H

#include <Eigen/Core>
#include <boost/scoped_array.hpp>
#include "diagonalizer_realhermitian.h"
#include "core.h"

class DiagonalizerRealHermitianLapack: public DiagonalizerRealHermitian
{
public:
	PROTEIN_CHAIN_API_DEBUG DiagonalizerRealHermitianLapack(unsigned int dim);
	PROTEIN_CHAIN_API_DEBUG unsigned int dim() const { return m_dim; }
	PROTEIN_CHAIN_API_DEBUG void diagonalize(const Eigen::MatrixXd& m, Eigen::MatrixXd& eigenvectors, Eigen::VectorXd& eigenvalues);
private:
	int m_dim;
	char m_jobz;
	int m_lwork;
	int m_liwork;
	int m_info;
	boost::scoped_array<double> m_work;
	boost::scoped_array<int> m_iwork;	
};

#endif // __DIAGONALIZER_REALHERMITIAN_LAPACK_H
