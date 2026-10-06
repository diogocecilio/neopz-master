/**
 * @file main.cpp
 * @brief Abaqus benchmark 1.15.2: consolidation of a triaxial specimen (Sects. 6.4 and 6.6).
 *
 * Usage: AbaqusTriaxialConsolidation [mp] [axi] [states] [3d] [novtk] (default: all the parts). The finite
 * element runs write the VTK file series of every increment in vtk/\<run name\>; "novtk" disables them.
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
