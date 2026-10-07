/**
 * @file main.cpp
 * @brief Embankment loading on a Modified Cam-Clay foundation: undrained loading and consolidation of a
 * 20 x 10 x 1 Hex20-Hex8 u-p slab in plane strain (Sect. 6.6, Figs. 12 to 14, Table 8 and the embankment
 * column of Table 10), with the elastic variant and the comparison of the tangent operators.
 *
 * Usage: EmbankmentConsolidation [novtk] [verbose] [model] [elastic] [tangents] [modes=D,sym,cont,DT,fd]
 *  - novtk: do not write the VTK file series of every converged state (vtk/embankment and
 *    vtk/embankment_elastic); the CSV files and the VTK files of the three states of Fig. 14 are always written;
 *  - verbose: one line per increment of the analyses;
 *  - model, elastic, tangents: run only these parts (default: all of them; see EmbankmentConsolidation::EPart);
 *  - modes=...: tangent operators of the comparison (default D,sym,cont,DT,fd).
 */
#include "EmbankmentConsolidation.h"

int main(int argc, char *argv[]) {
    EmbankmentConsolidation example;
    int parts = 0;
    for (int i = 1; i < argc; ++i) {
        const std::string arg(argv[i]);
        if (arg == "novtk") example.fWriteVTK = false;
        else if (arg == "verbose") example.fVerbose = 1;
        else if (arg == "model") parts |= EmbankmentConsolidation::EPartModel;
        else if (arg == "elastic") parts |= EmbankmentConsolidation::EPartElastic;
        else if (arg == "tangents") parts |= EmbankmentConsolidation::EPartTangents;
        else if (arg.rfind("modes=", 0) == 0) {
            example.fTangentModes.clear();
            std::stringstream list(arg.substr(6));
            std::string name;
            while (std::getline(list, name, ',')) {
                bool found = false;
                for (int m = 0; m <= TPZPlasticStepModifiedCamClay::EFiniteDifferenceTangent; ++m) {
                    const auto mode = static_cast<EmbankmentConsolidation::ETangentMode>(m);
                    if (name == TPZPlasticStepModifiedCamClay::TangentModeName(mode)) {
                        example.fTangentModes.push_back(mode);
                        found = true;
                    }
                }
                if (!found) {
                    std::cerr << "unknown tangent operator " << name << " (D, DT, sym, cont or fd)" << std::endl;
                    return 1;
                }
            }
        } else {
            std::cerr << "usage: " << argv[0] << " [novtk] [verbose] [model] [elastic] [tangents] [modes=D,sym,cont,DT,fd]"
                      << std::endl;
            return 1;
        }
    }
    if (parts) example.fParts = parts;
    example.RunAll();
    return 0;
}
