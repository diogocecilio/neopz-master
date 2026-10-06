/**
 * @file main.cpp
 * @brief Abaqus benchmark 1.15.2: consolidation of a triaxial specimen (Sects. 6.4 and 6.6).
 *
 * Usage: AbaqusTriaxialConsolidation [mp] [axi] [states] [3d] (default: all the parts).
 */
#include "AbaqusTriaxialConsolidation.h"

int main(int argc, char *argv[]) {
    std::vector<std::string> parts;
    for (int i = 1; i < argc; ++i) parts.push_back(argv[i]);
    if (parts.empty()) parts.push_back("all");
    AbaqusTriaxialConsolidation example;
    example.RunAll(parts);
    return 0;
}
