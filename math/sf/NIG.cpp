#include <cmath>
#include <boost/math/special_functions/bessel.hpp>
#include "NIG.h"

namespace rql
{
	namespace math
	{
		namespace sf
		{
			NIG::NIG(double alpha, double beta, double mu, double delta)
				: alpha_(alpha), beta_(beta), mu_(mu), delta_(delta), gamma_(gamma(alpha, beta))
			{
				mean_ = mean(beta, mu, delta, gamma_);
				variance_ = variance(alpha, delta, gamma_);
			}

			double NIG::gamma(double alpha, double beta)
			{
				return sqrt(alpha*alpha - beta*beta);
			}

			double NIG::mean(double beta, double mu, double delta, double gamma)
			{
				return mu + delta*beta/gamma;
			}

			double NIG::variance(double alpha, double delta, double gamma)
			{
				return delta*pow(alpha, 2)/pow(gamma, 3);
			}

			double NIG::pdf(double x) const
			{
				const double a = sqrt(delta_*delta_ + pow(x - mu_, 2));
				const double b = boost::math::cyl_bessel_k(1, alpha_ * a);	
				const double c = exp( delta_ * gamma_ + beta_ * (x - mu_) );
				return alpha_ * delta_ * b * c / a / 3.1415926535897932384626433832795028841971;
			}
		}
	}
}