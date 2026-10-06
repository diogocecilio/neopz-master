/**
 * @file main.cpp
 * @brief Embankment loading on a Modified Cam-Clay foundation: undrained loading and consolidation of a
 * plane strain u-p model (Sect. 6.5, Figs. 10 to 12 and Table 7), with the elastic and transposed-tangent
 * variants of the Python code.
 *
 * Usage: EmbankmentConsolidation [novtk]. By default the VTK file series of every converged state are
 * written in vtk/embankment and vtk/embankment_elastic; "novtk" writes only the CSV files and the VTK files
 * of the three states of Fig. 12.
 */
#include "EmbankmentConsolidation.h"

int main(int argc, char *argv[]) {
    EmbankmentConsolidation example;
    for (int i = 1; i < argc; ++i)
        if (std::string(argv[i]) == "novtk") example.fWriteVTK = false;
    example.RunAll();
    return 0;
}
