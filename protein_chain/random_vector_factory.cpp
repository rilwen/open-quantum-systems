#include "random_vector_factory.h"
#include "random_vector_gaussian.h"
#include "random_vector_levy_symmetric.h"
#include <boost/make_shared.hpp>
#include <cassert>

boost::shared_ptr<RandomVector> RandomVectorFactory::build_alpha_stable_impl(size_t len, double sigma, double alpha, const bool* antithetic)
{
	assert(alpha <= 2);
	assert(alpha > 0);
	if (alpha == 2) {
		if (antithetic) {
			return boost::make_shared<RandomVectorGaussian>(len, sigma, RandomVectorGaussian::DEFAULT_BATCH_SIZE, *antithetic);
		} else {
			return boost::make_shared<RandomVectorGaussian>(len, sigma);
		}
	} else {
		if (antithetic) {
			return boost::make_shared<RandomVectorLevySymmetric>(len, sigma, alpha, *antithetic);
		} else {
			return boost::make_shared<RandomVectorLevySymmetric>(len, sigma, alpha);
		}
	}
}
