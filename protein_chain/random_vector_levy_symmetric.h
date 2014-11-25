#ifndef __RANDOM_VECTOR_LEVY_SYMMETRIC_H
#define __RANDOM_VECTOR_LEVY_SYMMETRIC_H

#include "random_vector.h"
#include "core.h"
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_01.hpp>
#include <boost/random/variate_generator.hpp>

class RandomVectorLevySymmetric: public RandomVector
{
public:
	PROTEIN_CHAIN_API RandomVectorLevySymmetric(unsigned int dim, double sigma, double alpha, bool antithetic = false);
	PROTEIN_CHAIN_API_DEBUG virtual unsigned int dim() const { return m_dim; }
	PROTEIN_CHAIN_API_DEBUG unsigned int cutcount() const { return m_cutcount; };
	PROTEIN_CHAIN_API void draw(std::vector<double>& vec);
	PROTEIN_CHAIN_API_DEBUG void draw(Eigen::VectorXd& vec);
	PROTEIN_CHAIN_API_DEBUG void set_seed(unsigned int seed);
private:
	template <class V> void drawImpl(V& vec);
	double draw();
	void draw_batch();
private:
	static const unsigned int m_batchsize = 5000;
	const unsigned int m_dim;
	boost::mt19937 m_engine;
	boost::uniform_01<double> m_uniform_01;
	boost::variate_generator<boost::mt19937&, boost::uniform_01<double> > m_uniform_rng;
	double m_gamma;
	double m_alpha;
	std::vector<double> m_batch;
	std::vector<double>::const_iterator m_batch_iter;
	bool m_negate;
	unsigned int m_cutcount;
	bool m_antithetic;
};


#endif
