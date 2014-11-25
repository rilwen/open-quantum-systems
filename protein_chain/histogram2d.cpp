#include "histogram2d.h"
#include "interpolator2d.h"
#include <limits>
#include <stdexcept>
#include <cmath>

/*
static void find_boundaries(const std::vector<double>& v, double& min_v, double& max_v)
{
	min_v = std::numeric_limits<double>::infinity();
	max_v = -std::numeric_limits<double>::infinity();
	for (std::vector<double>::const_iterator it = v.begin(); it != v.end(); ++it) {
		min_v = std::min( *it, min_v );
		max_v = std::max( *it, max_v );
	}
}
*/

/*
Histogram2D::Histogram2D(const std::vector<double>& x, const std::vector<double>& y, unsigned int nx, unsigned int ny)
	: m_xpts(nx), m_ypts(ny), m_nbr_total(x.size()), m_hist(nx, ny)
{	
	if (x.size() != y.size()) {
		throw std::domain_error("x.size() != y.size()");
	}
	if (!x.size()) {
		throw std::domain_error("no data");
	}
	find_boundaries(x, m_min_x, m_max_x);
	find_boundaries(y, m_min_y, m_max_y);	
	init(nx, ny);
	fill(x, y);
}
*/

Histogram2D::Histogram2D(const std::vector<double>& x, const std::vector<double>& y, double minx, double maxx, unsigned int nx, double miny, double maxy, unsigned int ny)
	: m_xpts(nx), m_ypts(ny), m_nbr_total(x.size()), m_hist(nx, ny)
{	
	if (x.size() != y.size()) {
		throw std::domain_error("x.size() != y.size()");
	}
	if (!x.size()) {
		throw std::domain_error("no data");
	}
	if (minx >= maxx || miny >= maxy) {
		throw std::domain_error("wrong boundaries");
	}
	/*
	find_boundaries(x, m_min_x, m_max_x);
	find_boundaries(y, m_min_y, m_max_y);	
	m_min_x = std::max(m_min_x, minx);
	m_min_y = std::max(m_min_y, miny);
	m_max_x = std::min(m_max_x, maxx);
	m_max_y = std::min(m_max_y, maxy);
	*/
	m_min_x = minx;
	m_min_y = miny;
	m_max_x = maxx;
	m_max_y = maxy;
	init(nx, ny);
	fill(x, y);
}

void Histogram2D::init(unsigned int nx, unsigned int ny)
{
	m_hist.fill(0.0);
	if (nx == 0 || ny == 0) {
		throw std::domain_error("nx and ny must be > 0");
	}	
	m_dx = (m_max_x - m_min_x) / nx;
	m_dy = (m_max_y - m_min_y) / ny;
	assert( m_dx >= 0 );
	assert( m_dy >= 0 );
	if (m_dx == 0 || m_dy == 0) {
		throw std::runtime_error("point-like distribution");
	}
	for (unsigned int i = 0; i < m_xpts.size(); ++i) {
		m_xpts[i] = (i + 0.5) * m_dx + m_min_x;
	}
	for (unsigned int i = 0; i < m_ypts.size(); ++i) {
		m_ypts[i] = (i + 0.5) * m_dy + m_min_y;
	}
}

void Histogram2D::fill(const std::vector<double>& x, const std::vector<double>& y)
{
	const unsigned int nx = m_hist.rows();
	const unsigned int ny = m_hist.cols();
	std::vector<double>::const_iterator itx = x.begin();
	std::vector<double>::const_iterator ity = y.begin();
	assert( x.size() == y.size() );
	m_hist.fill(0.0);
	m_nbr_total = x.size();
	while (itx != x.end()) {
		assert( ity != y.end() );
		const double cx = *itx;
		const double cy = *ity;
		if (cx > m_min_x && cx < m_max_x && cy > m_min_y && cy < m_max_y) {
			unsigned int ix = std::min(static_cast<unsigned int>(floor((cx - m_min_x) / m_dx)), nx);
			unsigned int iy = std::min(static_cast<unsigned int>(floor((cy - m_min_y) / m_dy)), ny);
			++m_hist(ix, iy);
		}		
		++itx;
		++ity;
	}
	m_hist /= m_nbr_total;
}

/*
Interpolator2D Histogram2D::interpolator(double xShift, double xScale, double yShift, double yScale) const
{
	Interpolator2D rip(m_xpts, m_ypts, m_hist);	
	rip.shift()[0] = xShift;
	rip.shift()[1] = yShift;
	rip.scale()[0] = xScale;
	rip.scale()[1] = yScale;	
	const double factor = xScale*yScale;
	if (factor != 1) {
		rip *= factor / (m_dx*m_dy);
	} else {
		rip /= m_dx*m_dy;
	}
	return rip;
}
*/

/*
double Histogram2D::distance(const Histogram2D& left, double leftXShift, double leftXScale, double leftYShift, double leftYScale,
		const Histogram2D& right, double rightXShift, double rightXScale, double rightYShift, double rightYScale, unsigned int nx, unsigned int ny, double p, bool greedy)
{
	if (nx < 2 || ny < 2) {
		throw std::domain_error("nx and ny must be >= 2");
	}
	Interpolator2D leftInterp(left.interpolator(leftXShift, leftXScale, leftYShift, leftYScale));
	Interpolator2D rightInterp(right.interpolator(rightXShift, rightXScale, rightYShift, rightYScale));
	double min_x;
	double min_y;
	double max_x;
	double max_y;
	if (greedy) { // integrate over the union of the supports of the distributions
		min_x = std::min(leftInterp.minX(), rightInterp.minX());
		min_y = std::min(leftInterp.minY(), rightInterp.minY());
		max_x = std::max(leftInterp.maxX(), rightInterp.maxX());
		max_y = std::max(leftInterp.maxY(), rightInterp.maxY());
	} else { // integrate over the intersection of the supports of the distributions
		min_x = std::max(leftInterp.minX(), rightInterp.minX());
		min_y = std::max(leftInterp.minY(), rightInterp.minY());
		max_x = std::min(leftInterp.maxX(), rightInterp.maxX());
		max_y = std::min(leftInterp.maxY(), rightInterp.maxY());
	}
	return Interpolator2D::distance(leftInterp, rightInterp, min_x, max_x, nx, min_y, max_y, ny, p);
}

double Histogram2D::distance(const Histogram2D& left, const Histogram2D& right, unsigned int nx, unsigned int ny, double p, bool greedy)
{
	return distance(left, 0, 1, 0, 1, right, 0, 1, 0, 1, nx, ny, p, greedy);
}
*/
