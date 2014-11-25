#include "IntegrationDensities.h"
#include <cmath>
#include <stdexcept>

namespace rql {
	namespace integr {
		namespace IntegrationDensities {
			double Constant::operator()(double x0, double x1, unsigned int k) const
			{
				const double dx = x1 - x0;
				switch(k) {
				case 0:
					return m_a*dx;
				case 1:
					return m_a*dx*dx / 2;
				case 2:
					return m_a*dx*dx*dx / 6;
				case 3:
					return m_a*dx*dx*dx*dx / 24;
				default:
					throw std::domain_error("Not implemented");
				}
			}

			double Sine::operator()(double x0, double x1, unsigned int k) const
			{
				const double wx0 = m_omega*x0;
				const double wx1 = m_omega*x1;
				const double dwx = wx1 - wx0;
				const double c0 = cos(wx0);
				const double c1 = cos(wx1);
				const double s0 = sin(wx0);
				const double s1 = sin(wx1);
				switch(k) {
				case 0:
					return (c0 - c1) / m_omega;
				case 1:
					return (-dwx*c1 - s0 + s1) / m_omega / m_omega;
				case 2:
					return (-2*c0 - (dwx*dwx - 2)*c1 + 2*dwx*s1) / m_omega / m_omega / m_omega;
				case 3:
					return (-(dwx*dwx - 6)*dwx*c1 + 6*s0 + 3*(dwx*dwx - 2)*s1) / m_omega / m_omega / m_omega / m_omega;
				default:
					throw std::domain_error("Not implemented");
				}
			}

			double Cosine::operator()(double x0, double x1, unsigned int k) const
			{
				const double wx0 = m_omega*x0;
				const double wx1 = m_omega*x1;
				const double dwx = wx1 - wx0;
				const double c0 = cos(wx0);
				const double c1 = cos(wx1);
				const double s0 = sin(wx0);
				const double s1 = sin(wx1);
				switch(k) {
				case 0:
					return (s1 - s0) / m_omega;
				case 1:
					return (c1 - c0 + dwx*s1) / m_omega / m_omega;
				case 2:
					return (2*dwx*c1 + 2*s0 + (dwx*dwx - 2)*s1) / m_omega / m_omega / m_omega;
				case 3:
					return (6*c0 + 3*(dwx*dwx - 2)*c1 + dwx*s1*(dwx*dwx - 6)) / m_omega / m_omega / m_omega / m_omega;
				default:
					throw std::domain_error("Not implemented");
				}
			}
		}		
	}
}
