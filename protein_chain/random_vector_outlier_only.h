#ifndef __RANDOM_VECTOR_OUTLIER_ONLY_H
#define __RANDOM_VECTOR_OUTLIER_ONLY_H

#include "random_vector.h"
#include "core.h"
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_01.hpp>
#include <boost/random/variate_generator.hpp>

class RandomVectorOutlierOnly: public RandomVector
{
public:
	//! outlierProb in [0, 1]
	PROTEIN_CHAIN_API RandomVectorOutlierOnly(unsigned int dim, double outlierProb, double outlierSize);
	PROTEIN_CHAIN_API_DEBUG virtual unsigned int dim() const { return m_dim; }
	PROTEIN_CHAIN_API void draw(std::vector<double>& vec);
	PROTEIN_CHAIN_API_DEBUG void draw(Eigen::VectorXd& vec);
	PROTEIN_CHAIN_API_DEBUG void set_seed(unsigned int seed);
private:
	template <class V> void drawImpl(V& vec);
	const unsigned int m_dim;
	boost::mt19937 m_engine;
	boost::uniform_01<double> m_uniform_01;
	boost::variate_generator<boost::mt19937&, boost::uniform_01<double> > m_uniform_rng;
	double m_p;
	double m_p_h; // m_p / 2
	double m_outlier_size;
};

#endif // __RANDOM_VECTOR_OUTLIER_ONLY_H
