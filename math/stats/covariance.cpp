#include "covariance.h"

namespace rql { namespace math { namespace stats {

	Covariance::Covariance(unsigned int dim)
		: m_dim(dim), m_avg1(dim), m_avg2(boost::extents[dim][dim])
	{
	}

	Covariance::matrix_type Covariance::covariance() const
	{
		Eigen::MatrixXd cov(m_dim, m_dim);
		for (unsigned int i = 0; i < m_dim; ++i) {
			const double ai = m_avg1[i].value();
			cov(i, i) = m_avg2[i][i].value() - ai*ai;
			for (unsigned int j = i + 1; j < m_dim; ++j) {
				cov(j, i) = cov(i, j) = m_avg2[i][j].value() - ai*m_avg1[j].value();
			}
		}
		return cov;
	}

}}}