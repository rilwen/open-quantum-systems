#include "GSLFourierIntegral.h"

namespace rql {
	namespace math {
		namespace integr {
			GSLFourierIntegral::GSLFourierIntegral(double omega, double L, size_t n)
				: m_sin(omega, L, GSL_INTEG_SINE, n), m_cos(omega, L, GSL_INTEG_COSINE, n)
			{
			}
		}
	}
}
