#include "Jagged2DArrayTest.h"

typedef ::testing::Types<int, unsigned int, double> MyTypes;
INSTANTIATE_TYPED_TEST_CASE_P(My, Jagged2DArrayTest, MyTypes);
