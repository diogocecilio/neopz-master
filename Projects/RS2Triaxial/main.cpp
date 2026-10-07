/**
 * @file main.cpp
 * @brief Drained triaxial tests of the RS2 manual at a material point (Sect. 6.1, Fig. 5 and Table 2) and their
 * finite element check with one Hex20-Hex8 element (Fig. 4).
 */
#include "RS2Triaxial.h"

/** @brief Runs the example (the files are written to the current directory) */
int main() {
    RS2Triaxial example;
    example.RunAll();
    return 0;
}
