#ifndef __INTERPOLATOR_2D_H
#define __INTERPOLATOR_2D_H

#include <vector>
#include <Eigen/Core>
#include <boost/array.hpp>
#include "core.h"

class Interpolator2D
{
public:
	PROTEIN_CHAIN_API_DEBUG Interpolator2D(const std::vector<double>& x, const std::vector<double>& y, const Eigen::MatrixXd& z);
	// return 0 if outside bounds
	PROTEIN_CHAIN_API_DEBUG double operator()(double x, double y) const;
	PROTEIN_CHAIN_API_DEBUG boost::array<double,2>& shift() { return m_shift; }
	PROTEIN_CHAIN_API_DEBUG boost::array<double,2>& scale() { return m_scale; }
	PROTEIN_CHAIN_API_DEBUG const boost::array<double,2>& shift() const { return m_shift; }
	PROTEIN_CHAIN_API_DEBUG const boost::array<double,2>& scale() const { return m_scale; }
	PROTEIN_CHAIN_API_DEBUG double minX() const { return (m_x.front() - m_shift[0]) / m_scale[0]; }
	PROTEIN_CHAIN_API_DEBUG double maxX() const { return (m_x.back() - m_shift[0]) / m_scale[0]; }
	PROTEIN_CHAIN_API_DEBUG double minY() const { return (m_y.front() - m_shift[1]) / m_scale[1]; }
	PROTEIN_CHAIN_API_DEBUG double maxY() const { return (m_y.back() - m_shift[1]) / m_scale[1]; }
	PROTEIN_CHAIN_API_DEBUG Interpolator2D& operator*=(double s);
	PROTEIN_CHAIN_API_DEBUG Interpolator2D& operator/=(double s);
	PROTEIN_CHAIN_API_DEBUG static double distance(const Interpolator2D& left, const Interpolator2D& right, double minX, double maxX, unsigned int nx, double minY, double maxY, unsigned int ny, double p);
	PROTEIN_CHAIN_API_DEBUG double norm(unsigned int nx, unsigned int ny, double p) const;
private:
	//! Find such i that arr[i] <= val < arr[i+1], using binary search, assuming arr[arr.size()] == +infinity
	//! @return 0 <= i < arr.size()
	static unsigned int binary_search(const std::vector<double>& arr, double val);
private:
	const unsigned int m_dim_x;
	const unsigned int m_dim_y;
	std::vector<double> m_x;
	std::vector<double> m_y;
	Eigen::MatrixXd m_z;
	boost::array<double,2> m_shift;
	boost::array<double,2> m_scale;
	Eigen::MatrixXd m_slope_x;
	Eigen::MatrixXd m_slope_xy;
	Eigen::MatrixXd m_slope_y;
};

#endif // __INTERPOLATOR_2D_H
