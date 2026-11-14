#include "DrawdownTransient.h"

#include <algorithm>

DrawdownTransient::DrawdownTransient(TPZUnsatAnalysis &an, const DrawdownData &d, int idMacico,
                                     const std::set<int> &idsReservatorio)
    : fAn(an), fData(d)
{
    fFixed = an.AllDisplacements(TPZFMatrix<STATE>(an.NEq(), 1, 0.)); // esqueleto rigido: todo u = 0
    for (const auto &[no, eq] : an.AirEq()) fFixed[eq] = 0.;           // ar passivo
    for (const auto &[no, L] : an.BoundaryNodes(idMacico)) {           // lateral do macico abaixo do NA
        const REAL y = an.Coord(no, 1);
        if (y <= d.nivelMacico) fFixed[an.WaterEq().at(no)] = d.gammaw * (d.nivelMacico - y);
    }
    for (int id : idsReservatorio)
        for (const auto &[no, L] : an.BoundaryNodes(id)) fReservatorio.insert(no);
}

void DrawdownTransient::InitialState()
{
    TPZFMatrix<STATE> x(fAn.NEq(), 1, 0.);
    for (const auto &[no, eq] : fAn.WaterEq()) x(eq, 0) = fData.gammaw * (fData.nivelInicial - fAn.Coord(no, 1));
    fAn.X() = x;
    fAn.Restore();
    fAn.Accept(); // a memoria guarda p_w do passo anterior (o armazenamento parte dela)
    fFace.clear();
    fTime = 0.;
}

bool DrawdownTransient::AdvanceTo(REAL tHoras)
{
    const REAL alvo = tHoras * 3600., dt0 = fData.dt * 3600.;
    while (fTime < alvo - 1.e-6 * dt0) {
        const REAL dt = std::min(dt0, alvo - fTime);
        fTime += dt; // o passo usa o nivel do fim do passo
        if (!Step(dt)) {
            fTime -= dt;
            return false;
        }
    }
    return true;
}

bool DrawdownTransient::Step(REAL dt)
{
    const REAL tolP = 1.e-6, tolQ = 1.e-10;
    const REAL nivel = Level();
    // nos do reservatorio abaixo do nivel: hidrostatica; acima: candidatos a face de percolacao
    std::map<int64_t, REAL> submersos;
    std::set<int64_t> acima;
    for (int64_t no : fReservatorio) {
        const REAL y = fAn.Coord(no, 1);
        if (y <= nivel) submersos[fAn.WaterEq().at(no)] = fData.gammaw * (nivel - y);
        else acima.insert(no);
    }
    std::set<int64_t> face; // tentativa: a face anterior e os nos de cima com p_w > 0
    for (int64_t no : acima)
        if (fFace.count(no) || fAn.X().GetVal(fAn.WaterEq().at(no), 0) > tolP) face.insert(no);
    for (int it = 0; it < 40; it++) {
        std::map<int64_t, REAL> fixed = fFixed;
        for (const auto &[eq, v] : submersos) fixed[eq] = v;
        for (int64_t no : face) fixed[fAn.WaterEq().at(no)] = 0.;
        const TPZUnsatAnalysis::StepResult r = fAn.Solve(dt, fixed, {}, false);
        if (!r.conv) break;
        // confere a face: num no com p_w = 0 imposta, R_w > 0 seria agua entrando (sai da face); num no livre, p_w > 0
        // e agua empocando (entra na face)
        const TPZFMatrix<STATE> R = fAn.FullResidual();
        const TPZFMatrix<STATE> x = fAn.Mesh()->Solution();
        std::set<int64_t> nova = face;
        for (int64_t no : acima) {
            const int64_t eq = fAn.WaterEq().at(no);
            if (face.count(no)) {
                if (R.GetVal(eq, 0) / dt > tolQ) nova.erase(no);
            } else if (x.GetVal(eq, 0) > tolP) {
                nova.insert(no);
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

std::vector<REAL> DrawdownTransient::NodalWaterPressure() const
{
    std::vector<REAL> pw(fAn.Mesh()->Reference()->NNodes(), 0.);
    for (const auto &[no, eq] : fAn.WaterEq()) pw[no] = fAn.X().GetVal(eq, 0);
    return pw;
}
