#include "random_vector_levy_symmetric.h"
#include <stdexcept>
#include <cassert>
#include <cmath>

RandomVectorLevySymmetric::RandomVectorLevySymmetric(unsigned int dim, double sigma, double alpha, bool antithetic)
	: m_dim(dim), m_engine(17), m_uniform_01(), m_uniform_rng(m_engine, m_uniform_01)
	, m_gamma(alpha != 2 ? sigma : sigma/sqrt(2.0)) // !!!! so that sigma = 1 and alpha = 2 give us N(0,1) and for other alpha sigma is the scale parameter (consistent with other authors)
	, m_alpha(alpha), m_batch(m_batchsize * dim), m_batch_iter(m_batch.end()), m_negate(false), m_cutcount(0), m_antithetic(antithetic)
{	
	if (alpha <= 0 || alpha > 2) throw std::domain_error("COMPUTER SAYS NO");
}

template <class V> void RandomVectorLevySymmetric::drawImpl(V& vec)
{
	assert( vec.size() == m_dim );
	if (m_batch_iter == m_batch.end()) {
		draw_batch();
		m_batch_iter = m_batch.begin();
	}
	for (size_t i = 0; i < m_dim; ++i) {
		assert( m_batch_iter != m_batch.end() );
		vec[i] = *m_batch_iter*m_gamma;
		++m_batch_iter;
	}
}

void RandomVectorLevySymmetric::draw(std::vector<double>& vec)
{
	drawImpl(vec);
}

void RandomVectorLevySymmetric::draw(Eigen::VectorXd& vec)
{
	drawImpl(vec);
}

void RandomVectorLevySymmetric::draw_batch()
{
	static const double cutoff = 3000;
	if (m_negate) {
		for (std::vector<double>::iterator it = m_batch.begin(); it != m_batch.end(); ++it) {			
			*it *= -1;		
		}
		m_negate = false;
	} else {
		for (std::vector<double>::iterator it = m_batch.begin(); it != m_batch.end(); ++it) {
			*it = draw();
//			while (*it > cutoff || *it < -cutoff) {*it = draw(); m_cutcount++;}
		}
		if (m_antithetic) {
			m_negate = true;
		}
	}
}

static const double pi = 3.141592653589793;

double RandomVectorLevySymmetric::draw()
{
	if (m_alpha == 1) {
		// r = tan( pi/2 * (2*rand(sizeOut) - 1) ); 
		const double r = tan(0.5*pi*(2*m_uniform_rng() - 1));
		return r;
	} else {
		const double V = pi/2 * (2*m_uniform_rng() - 1);
		const double W = -log(m_uniform_rng());
		const double r = sin(m_alpha * V) / ( pow(cos(V),(1/m_alpha)) ) * pow( cos( V*(1-m_alpha) ) / W , (1-m_alpha)/m_alpha );

		return r;
	}
}

void RandomVectorLevySymmetric::set_seed(unsigned int seed)
{
	m_engine.seed(seed);
}

