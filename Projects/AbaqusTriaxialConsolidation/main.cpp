/**
 * @file main.cpp
 * @brief Abaqus benchmark 1.15.2: consolidation of a triaxial specimen with the Hex20-Hex8 model (Sects. 6.5
 * and 6.7 of the article).
 *
 * Usage: AbaqusTriaxialConsolidation [mesh] [mp] [fe] [tangents] [tolerance] [states] [softening] [novtk]
 * (default: all the parts).
 * The finite element runs of the parts fe, states and softening write the VTK file series of their increments in
 * vtk/\<run name\>; "novtk" disables them.
 */
#include "AbaqusTriaxialConsolidation.h"

int main(int argc, char *argv[]) {
    AbaqusTriaxialConsolidation example;
    std::vector<std::string> parts;
    for (int i = 1; i < argc; ++i) {
        if (std::string(argv[i]) == "novtk") example.fWriteVTK = false;
        else parts.push_back(argv[i]);
    }
    if (parts.empty()) parts.push_back("all");
    example.RunAll(parts);
    return 0;
}
