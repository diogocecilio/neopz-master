/**
 * @file main.cpp
 * @brief Terzaghi's consolidation of an elastic column with 1 x 1 x 10 Hex20-Hex8 elements (Sect. 6.4, Fig. 7,
 * Table 6).
 */
#include "TerzaghiConsolidation.h"

int main() {
    TerzaghiConsolidation example;
    example.RunAll();
    return 0;
}
