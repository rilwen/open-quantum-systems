#include "random_vector_threadsafe.h"
#include "random_vector.h"

RandomVectorThreadsafe::RandomVectorThreadsafe(boost::shared_ptr<RandomVector> rv)
	: m_rv(rv)
{
}

void RandomVectorThreadsafe::draw(std::vector<double>& epsilons)
{
	boost::lock_guard<boost::mutex> lock(m_mtx);
	m_rv->draw(epsilons);
}
