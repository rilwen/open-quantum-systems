#ifndef __HISTOGRAM_GENERATOR_H
#define __HISTOGRAM_GENERATOR_H

#include <boost/array.hpp>
#include <boost/multi_array.hpp>

// NOT FINISHED
template <unsigned int D> class HistogramGenerator
{
public:
	HistogramGenerator(const boost::array<double,D>& lb, const boost::array<double,D>& ub, const boost::array<unsigned int,D>& nbins);
	void update(const boost::array<double,D>& sample);
private:
	const boost::array<double,D> m_lb;
	const boost::array<double,D> m_ub;
	const boost::array<unsigned int,D> m_nbins;
	double m_total_cnt;
	boost::multi_array<double,D> m_cntrs;
};

template <unsigned int D> HistogramGenerator<D>::HistogramGenerator(const boost::array<double,D>& lb, const boost::array<double,D>& ub, const boost::array<unsigned int,D>& nbins)
	: m_lb(lb), m_ub(ub), m_nbins(nbins), m_total_cnt(0.0)
{
}

#endif // __HISTOGRAM_GENERATOR_H
