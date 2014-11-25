#include "random_vector_outlier_only.h"
#include <stdexcept>

RandomVectorOutlierOnly::RandomVectorOutlierOnly(unsigned int dim, double outlierProb, double outlierSize)
	: m_dim(dim), m_engine(17), m_uniform_01(), m_uniform_rng(m_engine, m_uniform_01), m_p(outlierProb), m_p_h(outlierProb/2), m_outlier_size(outlierSize)
{
	if (outlierProb < 0 || outlierProb > 1) {
		throw std::domain_error("Outlier probability outside [0, 1]");
	}
	if (outlierSize < 0) {
		throw std::domain_error("Outlier size negative");
	}
}

template <class V> void RandomVectorOutlierOnly::drawImpl(V& vec)
{
	for (unsigned int i = 0; i < m_dim; ++i) {
		const double x = m_uniform_rng();
		if (x < m_p_h) {
			vec[i] = m_outlier_size;
		} else if (x < m_p) {
			vec[i] = -m_outlier_size;
		} else {
			vec[i] = 0;
		}
	}
}

void RandomVectorOutlierOnly::draw(std::vector<double>& vec)
{
	drawImpl(vec);
}

void RandomVectorOutlierOnly::draw(Eigen::VectorXd& vec)
{
	drawImpl(vec);
}

void RandomVectorOutlierOnly::set_seed(unsigned int seed)
{
	m_engine.seed(seed);
}
