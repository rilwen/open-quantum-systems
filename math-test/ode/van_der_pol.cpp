#include "van_der_pol.h"
#include <gsl/gsl_errno.h>

VanDerPol::VanDerPol(double mu)
	: m_mu(mu)
{
}

int VanDerPol::operator()(double t, const double* y, double* f) const
{
	f[0] = y[1];
	f[1] =  -y[0] - m_mu*y[1]*(y[0]*y[0] - 1);
	return GSL_SUCCESS;
}

VanDerPol::Vec VanDerPol::operator()(double t, const VanDerPol::Vec& y) const
{
	VanDerPol::Vec result;
	(*this)(t, y, result);
	return result;
}

void VanDerPol::operator()(double t, const VanDerPol::Vec& y, VanDerPol::Vec& result) const
{	
	(*this)(t, y.data(), result.data());	
}

const double VanDerPol::defaultTime = 100;

const double VanDerPol::defaultResult = -1.75888;
