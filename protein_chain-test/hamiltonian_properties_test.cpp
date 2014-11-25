#include <gtest/gtest.h>
#include <protein_chain/hamiltonian_properties.h>

class HamiltonianPropertiesTest: public testing::Test
{
};

TEST_F(HamiltonianPropertiesTest,absorption_band_shift_to_interation_strength_ring_nn)
{
	ASSERT_EQ(0, HamiltonianProperties::absorption_band_shift_to_interation_strength_ring_nn(0, -1.3));
	ASSERT_EQ(0, HamiltonianProperties::absorption_band_shift_to_interation_strength_ring_nn(0, 1.3));
	ASSERT_EQ(0, HamiltonianProperties::absorption_band_shift_to_interation_strength_ring_nn(1, -1.3));
	ASSERT_EQ(0, HamiltonianProperties::absorption_band_shift_to_interation_strength_ring_nn(1, 1.3));
	ASSERT_NEAR(-1.3, HamiltonianProperties::absorption_band_shift_to_interation_strength_ring_nn(2, -1.3), 1E-14);
	ASSERT_NEAR(1.3, HamiltonianProperties::absorption_band_shift_to_interation_strength_ring_nn(2, 1.3), 1E-14);
	ASSERT_NEAR(-1.3, HamiltonianProperties::absorption_band_shift_to_interation_strength_ring_nn(3, -2.6), 1E-14);
	ASSERT_NEAR(1.3, HamiltonianProperties::absorption_band_shift_to_interation_strength_ring_nn(3, 2.6), 1E-14);
	ASSERT_NEAR(-1.3, HamiltonianProperties::absorption_band_shift_to_interation_strength_ring_nn(21, -2.6), 1E-14);
	ASSERT_NEAR(1.3, HamiltonianProperties::absorption_band_shift_to_interation_strength_ring_nn(21, 2.6), 1E-14);
}
