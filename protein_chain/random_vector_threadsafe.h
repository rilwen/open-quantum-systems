#ifndef __RANDOM_VECTOR_THREADSAFE_H
#define __RANDOM_VECTOR_THREADSAFE_H

#include <boost/shared_ptr.hpp>
#include <boost/thread.hpp>
#include "core.h"

class RandomVector;

class RandomVectorThreadsafe
{
public:
	PROTEIN_CHAIN_API RandomVectorThreadsafe(boost::shared_ptr<RandomVector> rv);
	PROTEIN_CHAIN_API void draw(std::vector<double>& epsilons);
private:
	boost::shared_ptr<RandomVector> m_rv;
	boost::mutex m_mtx;
};

#endif // __RANDOM_VECTOR_THREADSAFE_H
