#include "SlopeRain.h"

#include <algorithm>

namespace slope {

RainTransient::RainTransient(TPZUnsatAnalysis &an, const RainData &r, REAL gammaw) : fAn(an), fRain(r), fGammaW(gammaw)
{
    fFixed = an.AllDisplacements(TPZFMatrix<STATE>(an.NEq(), 1, 0.)); // esqueleto rigido: todo u = 0
    for (const auto &[no, eq] : an.AirEq()) fFixed[eq] = 0.;           // ar passivo
    for (int id : {bcEsquerda, bcDireita})                             // laterais abaixo do nivel: hidrostatica
        for (const auto &[no, L] : an.BoundaryNodes(id)) {
            const REAL nivel = r.WaterLevel(an.Coord(no, 0)), y = an.Coord(no, 1);
            if (y <= nivel) fFixed[an.WaterEq().at(no)] = gammaw * (nivel - y);
        }
    fSurface = an.BoundaryNodes(bcSuperficie); // chuva por metro de superficie (inclinada na face)
}

void RainTransient::InitialState()
{
    TPZFMatrix<STATE> x(fAn.NEq(), 1, 0.);
    for (const auto &[no, eq] : fAn.WaterEq())
        x(eq, 0) = std::max(fGammaW * (fRain.WaterLevel(fAn.Coord(no, 0)) - fAn.Coord(no, 1)), -fRain.suctionMax);
    fAn.X() = x;
    fAn.Restore();
    fAn.Accept(); // a memoria guarda p_w do passo anterior (o armazenamento do primeiro passo parte dela)
    fFace.clear();
    fTime = 0.;
    fSteps = fNewton = 0;
}

bool RainTransient::AdvanceTo(REAL t)
{
    const REAL alvo = t * 3600., dt0 = fRain.dt * 3600., fim = fRain.duration * 3600.;
    while (fTime < alvo - 1.e-6 * dt0) {
        const REAL dt = std::min(dt0, alvo - fTime);
        // chuva media no passo (um passo que cruza o fim da chuva recebe so a parte de antes)
        const REAL q = fRain.Intensity() * std::max(REAL(0), std::min(fTime + dt, fim) - fTime) / dt;
        if (!Step(dt, q)) return false;
        fTime += dt;
        fSteps++;
    }
    return true;
}

bool RainTransient::Step(REAL dt, REAL q)
{
    const REAL tolP = 1.e-6, tolQ = 1.e-10;
    std::set<int64_t> face = fFace; // tentativa: a face do passo anterior e os nos que ja tem p_w > 0
    for (const auto &[no, L] : fSurface)
        if (fAn.Nodal(fAn.X(), no, 1) > tolP) face.insert(no);
    for (int it = 0; it < 40; it++) {
        std::map<int64_t, REAL> fixed = fFixed, chuva;
        for (const auto &[no, L] : fSurface) {
            const int64_t eq = fAn.WaterEq().at(no);
            if (face.count(no)) fixed[eq] = 0.; // face de percolacao
            else chuva[eq] = q * L;             // vazao nodal da chuva (m3/s por metro de espessura)
        }
        const TPZUnsatAnalysis::StepResult r = fAn.Solve(dt, fixed, chuva, false);
        fNewton += r.iter;
        if (!r.conv) break;
        // confere a face: num no com p_w = 0 imposta, R_w / dt e a vazao que entra por ele
        const TPZFMatrix<STATE> R = fAn.FullResidual();
        const TPZFMatrix<STATE> x = fAn.Mesh()->Solution();
        std::set<int64_t> nova = face;
        for (const auto &[no, L] : fSurface) {
            const int64_t eq = fAn.WaterEq().at(no);
            if (face.count(no)) {
                if (R.GetVal(eq, 0) / (dt * L) > q + tolQ) nova.erase(no); // absorveria mais que a chuva: volta ao fluxo
            } else if (x.GetVal(eq, 0) > tolP) {
                nova.insert(no); // p_w > 0 sob a chuva: a agua empoca, o no vai para a face
            }
        }
        if (nova == face) {
            fAn.Accept();
            fFace = face;
            return true;
        }
        face = nova;
    }
    fAn.Restore();
    return false;
}

std::pair<REAL, REAL> RainTransient::PressureRange() const
{
    REAL pmin = 1.e300, pmax = -1.e300;
    for (const auto &[no, eq] : fAn.WaterEq()) {
        const REAL p = fAn.X().GetVal(eq, 0);
        pmin = std::min(pmin, p), pmax = std::max(pmax, p);
    }
    return {pmin, pmax};
}

} // namespace slope
