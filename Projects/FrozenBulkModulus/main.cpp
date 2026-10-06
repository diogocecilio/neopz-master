/**
 * @file main.cpp
 * @brief Exact integration of the porous law versus bulk modulus frozen at its trial value in triaxial tests
 * at a material point (Sect. 6.1, Table 3).
 */
#include "FrozenBulkModulus.h"

int main() {
    FrozenBulkModulus example;
    example.RunAll();
    return 0;
}
