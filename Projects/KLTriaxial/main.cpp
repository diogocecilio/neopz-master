/**
 * @file main.cpp
 * @brief Triaxial tests of example 1 of Krabbenhoft and Lyamin (2012) with one Hex20-Hex8 u-p element: drained and
 * undrained, NC and OCR = 10, against the closed forms and the material point (Sect. 7.1 of the article, v1.0).
 */
#include "KLTriaxial.h"

/**
 * @brief Runs the example (the files are written to the current directory); with the argument "tangents" it runs
 * only the comparison of the tangent operators in the drained tests (KLTriaxial::RunTangents)
 */
int main(int argc, char *argv[]) {
    KLTriaxial example;
    if (argc > 1 && std::string(argv[1]) == "tangents") {
        example.RunTangents();
        return 0;
    }
    example.RunAll();
    return 0;
}
