#ifndef __THERMAL_STATES_H
#define __THERMAL_STATES_H

#include<Eigen/Core>
#include "core.h"

namespace ThermalStates
{
	//! rho := exp(-H/kBT) / Trace[ exp(-H/kBT) ]; returns Trace[ exp(-H/kBT) ]
	PROTEIN_CHAIN_API std::complex<double> thermal_density_matrix(const Eigen::MatrixXcd& H, double kBT, Eigen::MatrixXcd& rho);

	//! Average mode occupation in a thermal state with temperature T (given by kB*T)
	PROTEIN_CHAIN_API double average_mode_occupation(double omega, double kBT);

	PROTEIN_CHAIN_API double average_mode_energy(double omega, double kBT);	
}

#endif // __THERMAL_STATES_H
