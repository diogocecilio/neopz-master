#include "RainTransient.h"

#include <algorithm>
#include <cmath>
#include <cstdio>

RainTransient::RainTransient(TPZUnsatAnalysis &an, const RainData &r, REAL gammaw, const std::set<int> &laterais,
                             const std::set<int> &superficie, const std::map<int, REAL> &cargas)
    : fAn(an), fRain(r), fGammaW(gammaw)
{
    fFixed = an.AllDisplacements(TPZFMatrix<STATE>(an.NEq(), 1, 0.)); // esqueleto rigido: todo u = 0
    for (const auto &[no, eq] : an.AirEq()) fFixed[eq] = 0.;           // ar passivo
    for (int id : laterais) {                                          // laterais abaixo do nivel: hidrostatica
        const auto h = cargas.find(id);
        for (const auto &[no, L] : an.BoundaryNodes(id)) {
            const REAL nivel = h != cargas.end() ? h->second : r.WaterLevel(an.Coord(no, 0)), y = an.Coord(no, 1);
            if (y <= nivel) fFixed[an.WaterEq().at(no)] = gammaw * (nivel - y);
        }
    }
    // chuva por metro de superficie (inclinada na face; ou da projecao horizontal, chuva vertical); um no entre dois
    // trechos soma as metades das arestas dos dois
    for (int id : superficie)
        for (const auto &[no, L] : an.BoundaryNodes(id, r.rainVertical)) fSurface[no] += L;
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
    const REAL alvo = t * 3600., dt0 = fRain.dt * 3600.;
    while (fTime < alvo - 1.e-6 * dt0) {
        const REAL dt = std::min(dt0, alvo - fTime);
        if (!Advance(dt)) return false;
    }
    return true;
}

bool RainTransient::Advance(REAL dt)
{
    const REAL fim = fRain.duration * 3600.;
    // chuva media no passo (um passo que cruza o fim da chuva recebe so a parte de antes)
    const REAL q = fRain.Intensity() * std::max(REAL(0), std::min(fTime + dt, fim) - fTime) / dt;
    if (Step(dt, q)) {
        fTime += dt;
        fSteps++;
        return true;
    }
    if (!(fRain.dtMin > 0.) || 0.5 * dt < fRain.dtMin * 3600. * (1. - 1.e-9)) return false;
    return Advance(0.5 * dt) && Advance(0.5 * dt);
}

bool RainTransient::SteadyState(REAL q, REAL dt0, REAL tol, REAL *mudanca)
{
    REAL dt = dt0 * 3600., dp = 1.e300;
    const REAL dtFim = 1.e4 * 3600., dtMax = 1.e8 * 3600.;
    fSteadySteps = 0;
    for (int k = 0; k < 200; k++) {
        const TPZFMatrix<STATE> antes = fAn.X();
        if (!Step(dt, q)) {
            dt *= 0.5;
            if (dt < 1.) break;
            continue;
        }
        fSteadySteps++;
        dp = 0.;
        for (const auto &[no, eq] : fAn.WaterEq()) dp = std::max(dp, std::abs(fAn.X().GetVal(eq, 0) - antes.GetVal(eq, 0)));
        if (dt >= dtFim && dp < tol) {
            fTime = 0.;
            if (mudanca) *mudanca = dp;
            return true;
        }
        dt = std::min(2. * dt, dtMax);
    }
    fTime = 0.;
    if (mudanca) *mudanca = dp;
    return false;
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
        const TPZUnsatAnalysis::StepResult r = fAn.Solve(dt, fixed, chuva, false, fRain.lineSearch);
        fNewton += r.iter;
        if (!r.conv) {
            char b[160];
            std::snprintf(b, sizeof b, "Newton (dt %.5f h, %d iteracoes, |r_p| %.2e, %zu nos na face)", dt / 3600., r.iter,
                          r.normp, face.size());
            fFalha = b;
            fAn.Restore();
            return false;
        }
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
    fFalha = "face de percolacao sem repetir em 40 tentativas (dt " + std::to_string(dt / 3600.) + " h, " +
             std::to_string(face.size()) + " nos)";
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

std::vector<REAL> RainTransient::NodalWaterPressure() const
{
    std::vector<REAL> pw(fAn.Mesh()->Reference()->NNodes(), 0.);
    for (const auto &[no, eq] : fAn.WaterEq()) pw[no] = fAn.X().GetVal(eq, 0);
    return pw;
}
