#include "scaling_params.h"
#include <stdexcept>

namespace ScalingParams
{
	double alpha_exponent(double alpha)
	{
		return 2*alpha / (1 + alpha);
	}

	void energy(double alpha, double& a, double& b, double& alphaExponent)
	{
		alphaExponent = alpha_exponent(alpha);
		if (alpha == 2) { // Gauss scaling params:
			a = 0;
			b = 1;
		} else if (alpha == 1) { // Lorentz scaling params
			a = 1;
			b = 4;
		} else if (alpha == 0.5) { // Levy scaling params
			a = -0.202815;
			b = 2.20413;
		} else {
			throw std::domain_error("No scaling parameters for this alpha");
		}
	}
}
