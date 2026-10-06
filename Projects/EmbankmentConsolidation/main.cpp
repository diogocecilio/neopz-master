/**
 * @file main.cpp
 * @brief Embankment loading on a Modified Cam-Clay foundation: undrained loading and consolidation of a
 * plane strain u-p model (Sect. 6.5, Figs. 10 to 12 and Table 7), with the elastic and transposed-tangent
 * variants of the Python code.
 */
#include "EmbankmentConsolidation.h"

int main() {
    EmbankmentConsolidation example;
    example.RunAll();
    return 0;
}
