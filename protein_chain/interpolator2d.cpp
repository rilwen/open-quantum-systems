#include "interpolator2d.h"
#include <stdexcept>
#include <cassert>

Interpolator2D::Interpolator2D(const std::vector<double>& x, const std::vector<double>& y, const Eigen::MatrixXd& z)
	: m_dim_x(x.size()), m_dim_y(y.size()), m_x(x), m_y(y), m_z(z)
{
	m_scale.fill(1);
	m_shift.fill(0);
	if (m_x.size() < 2 || m_y.size() < 2) {
		throw std::domain_error("x and y vector must have size() >= 2");
	}
	if (m_z.rows() != m_x.size() || m_z.cols() != m_y.size()) {
		throw std::domain_error("z matrix dimensions do not match x and y");
	}
	for (unsigned int i = 1; i < m_x.size(); ++i) {
		if (m_x[i] <= m_x[i - 1]) {
			throw std::domain_error("x values must be strictly increasing");
		}
	}
	for (unsigned int i = 1; i < m_y.size(); ++i) {
		if (m_y[i] <= m_y[i - 1]) {
			throw std::domain_error("y values must be strictly increasing");
		}
	}
	m_slope_x.resize(x.size(), y.size());
	m_slope_y.resize(x.size(), y.size());
	m_slope_xy.resize(x.size(), y.size());

	// calculate the slope of z along x and y
	for (unsigned int ix = 0; ix < m_dim_x; ++ix) {
		for (unsigned int iy = 0; iy < m_dim_y; ++iy) {
			const double zcurr = m_z(ix, iy);
			if (ix < m_dim_x - 1) {
				assert( m_x[ix + 1] - m_x[ix] > 0 );
				m_slope_x(ix, iy) = (m_z(ix+1, iy) - zcurr) / (m_x[ix + 1] - m_x[ix]);
			} else {
				m_slope_x(ix, iy) = 0;
			}
			if (iy < m_dim_y - 1) {
				assert( m_y[iy + 1] - m_y[iy] > 0 );
				m_slope_y(ix, iy) = (m_z(ix, iy + 1) - zcurr) / (m_y[iy + 1] - m_y[iy]);
			} else {
				m_slope_y(ix, iy) = 0;
			}
		}
	}

	// calculate the slope along y of the slope of z along x
	for (unsigned int ix = 0; ix < m_dim_x; ++ix) {
		for (unsigned int iy = 0; iy < m_dim_y; ++iy) {
			if (ix < m_dim_x - 1 && iy < m_dim_y - 1) {
				m_slope_xy(ix, iy) = (m_slope_x(ix, iy + 1) - m_slope_x(ix, iy)) / (m_y[iy + 1] - m_y[iy]);
			} else {
				m_slope_xy(ix, iy) = 0;
			}
		}
	}
}

double Interpolator2D::operator()(const double cx, const double cy) const
{
	// scale and shift back to original position
	// cx == (x - m_shift[0])/m_scale[0];
	// cy == (y - m_shift[1])/m_scale[1];
	const double x = m_scale[0] * cx + m_shift[0];
	const double y = m_scale[1] * cy + m_shift[1];
	if (x < m_x.front() || x > m_x.back() || y < m_y.front() || y > m_y.back()) {
		return 0;
	}
	const unsigned int ix = binary_search(m_x, x);
	const unsigned int iy = binary_search(m_y, y);
	assert( ix < m_dim_x );
	assert( iy < m_dim_y );
	const double dy = y - m_y[iy];
	return m_z(ix, iy) + m_slope_y(ix, iy) * dy + (m_slope_x(ix, iy) + m_slope_xy(ix, iy) * dy) * (x - m_x[ix]);
}

unsigned int Interpolator2D::binary_search(const std::vector<double>& arr, const double val)
{
	const unsigned int size = arr.size();
	assert( size > 1 );
	assert( val >= arr.front() );
	assert( val <= arr.back() );
	// logically, we assume that data[size] == +Infinity
	unsigned int l = 0; // highest i s.t. data[i] <= x
	unsigned int r = size; // lowest i s.t. data[i] > x
	unsigned int d = size;
	while (d > 1u)
	{
		const unsigned int m = l + d/2;
		if (val < arr[m])
			r = m;
		else
			l = m;
		d = r - l;
	}
	return l;
}

Interpolator2D& Interpolator2D::operator*=(double s)
{
	m_z *= s;
	m_slope_x *= s;
	m_slope_y *= s;
	m_slope_xy *= s;
	return *this;
}

Interpolator2D& Interpolator2D::operator/=(double s)
{
	m_z /= s;
	m_slope_x /= s;
	m_slope_y /= s;
	m_slope_xy /= s;
	return *this;
}

double Interpolator2D::distance(const Interpolator2D& left, const Interpolator2D& right, double minX, double maxX, unsigned int nx, double minY, double maxY, unsigned int ny, double p)
{
	const double dx = (maxX - minX) / (nx-1);
	const double dy = (maxY - minY) / (ny-1);
	double sum = 0;
	for (unsigned int i = 0; i < nx; ++i) {
		const double x = minX + i * dx;
		for (unsigned int j = 0; j < ny; ++j) {
			const double y = minY + j * dy;
			const double lz = left(x, y);
			const double rz = right(x, y);
			const double diff = lz - rz;
			sum += pow(std::abs(diff), p)*dx*dy;
		}
	}
	return pow(sum, 1/p);
}

double Interpolator2D::norm(unsigned int nx, unsigned int ny, double p) const
{
	Interpolator2D clone(*this);
	clone *= 0;
	return distance(*this, clone, minX(), maxX(), nx, minY(), maxY(), ny, p);
}