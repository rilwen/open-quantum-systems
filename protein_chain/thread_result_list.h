#ifndef __THREAD_RESULT_LIST_H
#define __THREAD_RESULT_LIST_H

#include <list>
#include <utility>
#include <boost/shared_ptr.hpp>
#include <boost/thread/mutex.hpp>
#include <boost/thread/locks.hpp>

template <class ResultType> class ThreadResultList
{
public:	
	virtual void push_back(unsigned int thread_index, const ResultType& result) = 0;
};

template <class ResultType> class ThreadResultListImpl: public ThreadResultList<ResultType>
{
public:	
	typedef std::pair<unsigned int, ResultType> item_type;
	typedef typename std::list<item_type>::const_iterator const_iterator;
public:
	//! Synchronized.
	void push_back(unsigned int thread_index, const ResultType& result);

	//! Not synchronized. Use after all threads which write to the list have finished.
	const_iterator begin() const { return m_data.begin(); }

	//! Not synchronized. Use after all threads which write to the list have finished.
	const_iterator end() const { return m_data.end(); }
private:	
	std::list<item_type> m_data;
	boost::mutex m_mtx;
};

template <class RT> void ThreadResultListImpl<RT>::push_back(unsigned int thread_index, const RT& result)
{
	boost::lock_guard<boost::mutex> lock(m_mtx);
	m_data.push_back(item_type(thread_index, result));
}

#endif // __THREAD_RESULT_LIST_H
