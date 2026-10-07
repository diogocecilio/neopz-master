/**
 * @file main.cpp
 * @brief Drained and undrained triaxial tests of FLAC3D with one Hex20-Hex8 u-p element (Sect. 6.3, Fig. 6,
 * Table 5) and the global iterations with five tangent operators (Sect. 6.7, Table 10).
 */
#include "FLAC3DTriaxial.h"

int main() {
    FLAC3DTriaxial example;
    example.RunAll();
    return 0;
}
