#ifndef __DOUBLE_EXPONENTIAL_QUADRATURE_H
#define __DOUBLE_EXPONENTIAL_QUADRATURE_H

#include <vector>
#include "core.h"

class DoubleExponentialQuadrature 
{
public:
	PROTEIN_CHAIN_API static void build(double a, double b, unsigned int n, std::vector<double>& samplingpnts, std::vector<double>& weights);
};

#endif
