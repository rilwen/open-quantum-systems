#ifndef __HISTOGRAM_2D_H
#define __HISTOGRAM_2D_H

#include <vector>
#include <Eigen/Core>
#include "core.h"

class Interpolator2D;

class Histogram2D
{
public:
	//Histogram2D(const std::vector<double>& x, const std::vector<double>& y, unsigned int nx, unsigned int ny);
	PROTEIN_CHAIN_API Histogram2D(const std::vector<double>& x, const std::vector<double>& y, double minx, double maxx, unsigned int nx, double miny, double maxy, unsigned int ny);
	PROTEIN_CHAIN_API const Eigen::MatrixXd& histogramMatrix() const { return m_hist; }
	PROTEIN_CHAIN_API_DEBUG unsigned int totalNumber() const { return m_nbr_total; }
	PROTEIN_CHAIN_API_DEBUG const std::vector<double>& xPts() const { return m_xpts; }
	PROTEIN_CHAIN_API_DEBUG const std::vector<double>& yPts() const { return m_ypts; }
	PROTEIN_CHAIN_API double dx() const { return m_dx; }
	PROTEIN_CHAIN_API double dy() const { return m_dy; }
	//Interpolator2D interpolator(double xShift = 0, double xScale = 1, double yShift = 0, double yScale = 1) const;
	//! p-norm
	/*
	static double distance(const Histogram2D& left, double leftXShift, double leftXScale, double leftYShift, double leftYScale,
		const Histogram2D& right, double rightXShift, double rightXScale, double rightYShift, double rightYScale, unsigned int nx, unsigned int ny, double p, bool greedy);		
	static double distance(const Histogram2D& left, const Histogram2D& right, unsigned int nx, unsigned int ny, double p, bool greedy);
	*/
	PROTEIN_CHAIN_API double norm() const { return m_hist.array().abs().sum(); }
	PROTEIN_CHAIN_API double distance(const Histogram2D& other) const { return (m_hist - other.m_hist).array().abs().sum(); }
	PROTEIN_CHAIN_API double minX() const { return m_min_x; }
	PROTEIN_CHAIN_API double minY() const { return m_min_y; }
	PROTEIN_CHAIN_API double maxX() const { return m_max_x; }
	PROTEIN_CHAIN_API double maxY() const { return m_max_y; }
	PROTEIN_CHAIN_API_DEBUG void fill(const std::vector<double>& x, const std::vector<double>& y);
	PROTEIN_CHAIN_API_DEBUG Histogram2D& operator*=(double scale) { m_hist *= scale; return *this; }
private:
	void init(unsigned int nx, unsigned int ny);		
private:
	double m_min_x;
	double m_max_x;
	double m_dx;
	double m_min_y;
	double m_max_y;
	double m_dy;
	std::vector<double> m_xpts;
	std::vector<double> m_ypts;
	unsigned int m_nbr_total;
	Eigen::MatrixXd m_hist;
};

#endif //__HISTOGRAM_2D_H
