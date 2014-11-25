#include "double_exponential_quadrature.h"
#include <cmath>

static const double A = 3.5;
static const double pi = 3.141592653589793;

void DoubleExponentialQuadrature::build(double a, double b, unsigned int n, std::vector<double>& samplingpnts, std::vector<double>& weights)
{
	samplingpnts.resize(n+1);
	weights.resize(n+1);
	const double h = 2*A/n;

	for (unsigned int k = 0; k <= n; ++k) {
		const double t = -A + k*h;
		//Quadrature over [-1,1] interval
		samplingpnts[k] = tanh(pi*sinh(t)/2);
		weights[k] = h*pi*cosh(t)/(1 + cosh(pi*sinh(t)));
		if (k==0 || k==n) weights[k] *= 0.5;
		//Rescale on the interval [a,b]
		weights[k] *= (b - a)/2;
		samplingpnts[k] = a + (samplingpnts[k] + 1)/2*(b - a);
	}
	
}
