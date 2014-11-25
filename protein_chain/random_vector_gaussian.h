#ifndef __RANDOM_VECTOR_GAUSSIAN_H
#define __RANDOM_VECTOR_GAUSSIAN_H

#include "random_vector.h"
#include "core.h"
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_01.hpp>
#include <boost/random/variate_generator.hpp>


class RandomVectorGaussian: public RandomVector
{
public:
	static const size_t DEFAULT_BATCH_SIZE;
	PROTEIN_CHAIN_API RandomVectorGaussian(unsigned int dim, double sigma = 1, size_t batch_size = DEFAULT_BATCH_SIZE, bool antithetic_sampling = true);
	PROTEIN_CHAIN_API_DEBUG unsigned int dim() const { return m_dim; }
	PROTEIN_CHAIN_API void draw(std::vector<double>& vec);
	PROTEIN_CHAIN_API_DEBUG void draw(Eigen::VectorXd& vec);
	PROTEIN_CHAIN_API_DEBUG void set_seed(unsigned int seed);
private:
	static double normsinv(double p);
	void draw_batch();
	void decorrelate();
	template <class V> void drawImpl(V& vec);
private:
	const unsigned int m_dim;
	boost::mt19937 m_engine;
	boost::uniform_01<double> m_uniform_01;
	boost::variate_generator<boost::mt19937&, boost::uniform_01<double> > m_uniform_rng;
	double m_sigma;
	std::vector<double> m_batch;
	std::vector<double>::const_iterator m_batch_iter;
	bool m_negate;
	bool m_antithetic_sampling;
};

#endif // __RANDOM_VECTOR_GAUSSIAN_H
