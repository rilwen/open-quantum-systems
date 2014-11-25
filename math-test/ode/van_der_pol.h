#ifndef __MATH_TEST_VAN_DER_POL_H
#define __MATH_TEST_VAN_DER_POL_H

#include <Eigen/Core>

class VanDerPol
{
public:
	typedef Eigen::Vector2d Vec;
	VanDerPol(double mu = 10.0);
	size_t dim() const { return 2; }
	int operator()(double t, const double* y, double* f) const;
	void operator()(double t, const Vec& y, Vec& f) const;
	Vec operator()(double t, const Vec& y) const;
	template <class V> static void initialConditions(V& y);
	static const double defaultTime;
	static const double defaultResult;
private:
	double m_mu;
};

template <class V> void VanDerPol::initialConditions(V& y)
{
	y[0] = 1;
	y[1] = 0;
}

#endif // __MATH_TEST_VAN_DER_POL_H
