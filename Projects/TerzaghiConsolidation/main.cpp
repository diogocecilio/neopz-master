/**
 * @file main.cpp
 * @brief Terzaghi's consolidation of an elastic column with 1 x 1 x 10 Hex20-Hex8 elements (Sect. 6.4, Fig. 8,
 * Table 6).
 */
#include "TerzaghiConsolidation.h"

/** @brief Runs the example (the files are written to the current directory) */
int main() {
    TerzaghiConsolidation example;
    example.RunAll();
    return 0;
}
