#ifndef __PROTEIN_CHAIN_RANDOM_VECTOR_FACTORY_H
#define __PROTEIN_CHAIN_RANDOM_VECTOR_FACTORY_H

#include <boost/shared_ptr.hpp>
#include "core.h"

class RandomVector;

class RandomVectorFactory
{
public:
	//! For backward compatibility: use default setting for antithetic sampling (TRUE for Gaussian, FALSE for Levy)
	PROTEIN_CHAIN_API static boost::shared_ptr<RandomVector> build_alpha_stable(size_t len, double sigma, double alpha)
	{
		return build_alpha_stable_impl(len, sigma, alpha, 0);
	}

	PROTEIN_CHAIN_API static boost::shared_ptr<RandomVector> build_alpha_stable(size_t len, double sigma, double alpha, bool antithetic)
	{
		return build_alpha_stable_impl(len, sigma, alpha, &antithetic);
	}
private:
	static boost::shared_ptr<RandomVector> build_alpha_stable_impl(size_t len, double sigma, double alpha, const bool* antithetic);
};


#endif // __PROTEIN_CHAIN_RANDOM_VECTOR_FACTORY_H
