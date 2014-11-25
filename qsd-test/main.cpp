#include <stdlib.h>
#include <gtest/gtest.h>
#include <iostream>

/*
 *
 */
int main(int argc, char** argv)
{
	try {
		testing::InitGoogleTest(&argc, argv);
		return RUN_ALL_TESTS();
	} catch (std::exception& e) {
		std::cerr << "Caught exception: " << e.what() << std::endl;
		return -1;
	}
}

