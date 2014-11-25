#include "hamiltonian_properties.h"
#include "hamiltonian_factory.h"
#include <vector>

namespace HamiltonianProperties {
	double absorption_band_shift_to_interation_strength_ring_nn(size_t dimension, double C)
	{
		if (dimension <= 1) return 0;
		if (dimension % 2 == 1)	return C/2;
		if (C==0) return 0;
		Eigen::MatrixXd h(dimension,dimension);
		if (C > 0) {
			HamiltonianFactory::buildRing(std::vector<double>(dimension, 0), 1, h);
			return C / highest_eigenenergy(h);
		} else {
			HamiltonianFactory::buildRing(std::vector<double>(dimension, 0), -1, h);
			return -C / lowest_eigenenergy(h);
		}
	}

	double absorption_band_shift_to_interation_strength_wheel_nn(size_t dimension, double C)
	{
		if (dimension <= 1)	return 0;
		if (C==0) return 0;
		Eigen::MatrixXd h(dimension,dimension);
		if (C > 0) {
			HamiltonianFactory::buildWheel(std::vector<double>(dimension, 0), 1, h);
			return C / highest_eigenenergy(h);
		} else {
			HamiltonianFactory::buildWheel(std::vector<double>(dimension, 0), -1, h);
			return -C / lowest_eigenenergy(h);
		}
	}
}
