#ifndef __SCALING_PARAMS_H
#define __SCALING_PARAMS_H

#include "core.h"

namespace ScalingParams
{
	PROTEIN_CHAIN_API double alpha_exponent(double alpha);

	//! Out of date!!!
	PROTEIN_CHAIN_API void energy(double alpha, double& a, double& b, double& alphaExponent);	
};

#endif // __SCALING_PARAMS_H
