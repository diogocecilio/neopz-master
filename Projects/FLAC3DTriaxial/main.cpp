/**
 * @file main.cpp
 * @brief Drained and undrained triaxial tests of FLAC3D with one Hex20-Hex8 u-p element (Sect. 6.3, Figs. 4 and 7,
 * Table 5) and the global iterations with five tangent operators (Sect. 6.7, Table 10).
 */
#include "FLAC3DTriaxial.h"

/** @brief Runs the example (the files are written to the current directory) */
int main() {
    FLAC3DTriaxial example;
    example.RunAll();
    return 0;
}
