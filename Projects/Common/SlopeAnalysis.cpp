#include "SlopeAnalysis.h"

#include "TPZSpStructMatrix.h"
#include "TPZUnsatAnalysis.h"
#include "pzstepsolver.h"
#include "Plasticity/TPZPoroElastoPlasticUPAnalysis.h"

namespace slope {

const char *Name(MaterialKind k) { return k == MaterialKind::UP ? "u-p" : "unsat"; }

const char *Name(WaterModel w)
{
    switch (w) {
    case WaterModel::Dry: return "seco";
    case WaterModel::Drawdown: return "rebaixamento, pressoes do MEF";
    case WaterModel::DrawdownForce: return "rebaixamento, forca -grad u do MEF";
    case WaterModel::ForceField: return "forca de percolacao dada";
    case WaterModel::Rain: return "chuva";
    }
    return "";
}

SlopeAnalysis::SlopeAnalysis(const Data &d, MaterialKind kind, WaterModel water, FactorMethod factor)
    : fData(d), fKind(kind), fWater(water), fFactor(factor)
{
    if (water == WaterModel::Rain) {
        std::cerr << "SlopeAnalysis: a chuva (transiente do material nao saturado) ainda nao esta implementada\n";
        DebugStop();
    }
    if (water == WaterModel::Rain && kind != MaterialKind::Unsat) DebugStop();
}

namespace {

void PrintLevelHeader(const std::string &header, int64_t neq, int64_t nfree)
{
    std::cout << header << " | " << neq << " equacoes (" << nfree << " livres)" << std::endl;
}

void PrintSeepage(const std::vector<REAL> &pw)
{
    Print("    percolacao (esqueleto rigido): p_w de %.1f a %.1f kPa", *std::min_element(pw.begin(), pw.end()),
          *std::max_element(pw.begin(), pw.end()));
}

// Propriedades de cada ponto de integracao (campos aleatorios)
template <class TMat>
void ApplyProperties(TPZMultiphysicsCompMesh *cmesh, TMat *mat, const SlopeAnalysis::PointProperties &f)
{
    if (!f) return;
    ForEachPoint(cmesh, [&](TPZCompEl *, const TPZVec<REAL> &x, int64_t im) { f(x, mat->MemItem(im)); });
}

// Forca de corpo em cada ponto de integracao: f(elemento, x)
template <class TMat>
void SetPointForces(TPZMultiphysicsCompMesh *cmesh, TMat *mat,
                    const std::function<std::array<REAL, 2>(TPZCompEl *, const TPZVec<REAL> &)> &f)
{
    auto forces = std::make_shared<std::vector<std::array<REAL, 2>>>(mat->GetMemory()->NElements(),
                                                                      std::array<REAL, 2>{0., 0.});
    ForEachPoint(cmesh, [&](TPZCompEl *cel, const TPZVec<REAL> &x, int64_t im) { (*forces)[im] = f(cel, x); });
    mat->fPointForce = forces;
}

// -grad u do campo P1 de p_w (por no) no triangulo do elemento, y para cima: (-dp/dx, -dp/dy - gamma_w)
std::array<REAL, 2> ForceFromNodal(TPZGeoEl *gel, const std::vector<REAL> &pw, REAL gammaw)
{
    REAL x[3], y[3], p[3];
    for (int k = 0; k < 3; k++) {
        const int64_t no = gel->NodeIndex(k);
        x[k] = gel->NodePtr(k)->Coord(0), y[k] = gel->NodePtr(k)->Coord(1), p[k] = pw[no];
    }
    const REAL dx2 = x[1] - x[0], dy2 = y[1] - y[0], dx3 = x[2] - x[0], dy3 = y[2] - y[0], det = dx2 * dy3 - dx3 * dy2;
    const REAL a = ((p[1] - p[0]) * dy3 - (p[2] - p[0]) * dy2) / det, b = (dx2 * (p[2] - p[0]) - dx3 * (p[1] - p[0])) / det;
    return {-a, -b - gammaw};
}

template <class TMat>
void WritePointForces(TPZMultiphysicsCompMesh *cmesh, TMat *mat, const std::string &file)
{
    if (!mat->fPointForce) return;
    std::ofstream out(file);
    out << "x,y,fx,fy\n";
    ForEachPoint(cmesh, [&](TPZCompEl *, const TPZVec<REAL> &x, int64_t im) {
        const auto &f = (*mat->fPointForce)[im];
        out << x[0] << "," << x[1] << "," << f[0] << "," << f[1] << "\n";
    });
}

} // namespace

SlopeAnalysis::LevelResult SlopeAnalysis::SolveLevelUP(TPZGeoMesh *gmesh, const std::vector<int> &ordens,
                                                       const std::string &header, const std::string &base)
{
    const Data &d = fData;
    const REAL peso = UsesForce() ? d.gamma - d.gammaw : d.gamma; // com a forca de percolacao: o peso submerso
    TPZMultiphysicsCompMesh *cmesh = CreateUPCompMesh(gmesh, 2, d, ordens);
    auto *mat = dynamic_cast<MatUP *>(cmesh->FindMaterial(matSolo));
    if (!mat) DebugStop();
    ApplyProperties(cmesh, mat, fProps);
    LevelResult res;
    {
        TPZPoroElastoPlasticUPAnalysis an(cmesh, mat);
        TPZSpStructMatrix<STATE> strmat(cmesh); // PARDISO, nao simetrica
        an.SetStructuralMatrix(strmat);
        TPZStepSolver<STATE> solver;
        solver.SetDirect(ELU);
        an.SetSolver(solver);
        an.IdentifyEquations(); // equacoes de pressao e com Dirichlet
        res.neq = cmesh->NEquations();
        std::vector<int64_t> livres; // do fator: deslocamentos sem apoio
        for (int64_t eq = 0; eq < res.neq; eq++)
            if (!an.PressureEquations()[eq] && !an.ConstrainedEquations()[eq]) livres.push_back(eq);
        res.nfree = livres.size();
        PrintLevelHeader(header, res.neq, res.nfree);

        TPZFMatrix<STATE> x0(res.neq, 1, 0.);
        if (fWater == WaterModel::Drawdown || fWater == WaterModel::DrawdownForce) {
            std::set<int64_t> contorno; // percolacao: u = 0 e p do contorno fixos
            for (const auto &[no, v] : BoundaryPressures(gmesh, d)) {
                const int64_t eq = an.NodeEquation(no, 1, 0);
                x0(eq, 0) = v;
                contorno.insert(eq);
            }
            std::vector<int64_t> livresP;
            for (int64_t eq = 0; eq < res.neq; eq++)
                if (an.PressureEquations()[eq] && !contorno.count(eq)) livresP.push_back(eq);
            mat->SetTimeStep(1.); // liga o fluxo; sem armazenamento o passo ja e o permanente
            if (!SolveFree(an, livresP, x0)) DebugStop();
            mat->SetTimeStep(0.);
            for (int64_t no = 0; no < gmesh->NNodes(); no++) res.pw.push_back(x0(an.NodeEquation(no, 1, 0), 0));
            PrintSeepage(res.pw);
            if (fWater == WaterModel::DrawdownForce) { // a percolacao vira forca; pressoes nulas
                SetPointForces(cmesh, mat, [&](TPZCompEl *cel, const TPZVec<REAL> &) {
                    return ForceFromNodal(cel->Reference(), res.pw, d.gammaw);
                });
                x0.Zero();
            }
        } else if (fWater == WaterModel::ForceField) {
            if (!fForce) DebugStop();
            SetPointForces(cmesh, mat, [&](TPZCompEl *, const TPZVec<REAL> &x) { return fForce(x); });
        }
        PathProblem p;
        p.an = &an;
        p.x0 = x0;
        p.scaleFixed = true; // pressoes fixas = lambda p
        p.free = livres;
        p.setLambda = [mat, peso](REAL lam) {
            mat->SetBodyForce(TPZManVector<REAL, 3>({0., -lam * peso, 0.}));
            mat->fForceFactor = lam;
        };
        p.setReduction = [mat](REAL F) { mat->PlasticModel().SetStrengthReductionFactor(F); };
        p.accept = [&an]() { QuietCout q; an.AcceptSolution(); };
        SetFreeEquations(an, p.free);
        const FactorResult f = ComputeFactor(p, res.neq, fFactor);
        res.value = f.value, res.steps = f.steps, res.end = f.end;
    }
    res.fmc = MaxYield(mat, d);
    res.ind = PlasticIndicator(cmesh, mat);
    if (!base.empty()) {
        WritePlasticPoints(cmesh, mat, base + "_pontos.csv");
        WritePointForces(cmesh, mat, base + "_forca.csv");
    }
    EraseMesh(cmesh);
    return res;
}

SlopeAnalysis::LevelResult SlopeAnalysis::SolveLevelUnsat(TPZGeoMesh *gmesh, const std::vector<int> &ordens,
                                                          const std::string &header, const std::string &base)
{
    const Data &d = fData;
    const REAL peso = UsesForce() ? d.gamma - d.gammaw : d.gamma;
    TPZMultiphysicsCompMesh *cmesh = CreateUnsatCompMesh(gmesh, 2, d, peso, ordens);
    auto *mat = dynamic_cast<MatUnsat *>(cmesh->FindMaterial(matSolo));
    if (!mat) DebugStop();
    ApplyProperties(cmesh, mat, fProps);
    LevelResult res;
    {
        // analise antes da matriz: o construtor renumera as equacoes
        TPZUnsatAnalysis an(cmesh, {{bcBase, {true, true}}, {bcDireita, {true, false}}, {bcEsquerda, {true, false}}});
        TPZSpStructMatrix<STATE> strmat(cmesh);
        an.SetStructuralMatrix(strmat);
        TPZStepSolver<STATE> solver;
        solver.SetDirect(ELU);
        an.SetSolver(solver);
        res.neq = an.NEq();
        std::vector<int64_t> livres; // do fator: deslocamentos sem apoio (pressoes congeladas)
        const auto fixed = an.FrozenFixed(an.X());
        for (int64_t eq = 0; eq < res.neq; eq++)
            if (!fixed.count(eq)) livres.push_back(eq);
        res.nfree = livres.size();
        PrintLevelHeader(header, res.neq, res.nfree);

        TPZFMatrix<STATE> x0(res.neq, 1, 0.);
        if (fWater == WaterModel::Drawdown || fWater == WaterModel::DrawdownForce) {
            auto fixas = an.AllDisplacements(an.X()); // percolacao: u = 0, p_a = 0 e p_w do contorno fixos
            for (const auto &[no, eq] : an.AirEq()) fixas[eq] = 0.;
            for (const auto &[no, v] : BoundaryPressures(gmesh, d)) fixas[an.WaterEq().at(no)] = v;
            {
                QuietCout q;
                if (!an.Solve(1., fixas, {}, false).conv) DebugStop(); // sem armazenamento: o permanente
                an.Accept();
            }
            x0 = an.X();
            an.ResetMemory(); // o fator parte do estado sem tensoes (as propriedades dos pontos ficam)
            for (int64_t no = 0; no < gmesh->NNodes(); no++) res.pw.push_back(an.Nodal(x0, no, 1));
            PrintSeepage(res.pw);
            if (fWater == WaterModel::DrawdownForce) {
                SetPointForces(cmesh, mat, [&](TPZCompEl *cel, const TPZVec<REAL> &) {
                    return ForceFromNodal(cel->Reference(), res.pw, d.gammaw);
                });
                x0.Zero();
            }
        } else if (fWater == WaterModel::ForceField) {
            if (!fForce) DebugStop();
            SetPointForces(cmesh, mat, [&](TPZCompEl *, const TPZVec<REAL> &x) { return fForce(x); });
        }
        an.SetGravityScalesSuction(true); // lambda escala o peso e toda a pressao (o Gamma do artigo)
        for (auto *m : an.Materials()) {
            m->SetTimeStep(0.);
            m->SetFrozenPressures(true); // so o equilibrio
        }
        PathProblem p;
        p.an = &an;
        p.x0 = x0;
        p.scaleFixed = false; // o material escala as pressoes congeladas
        p.free = livres;
        p.setLambda = [&an, mat](REAL lam) {
            an.SetGravityFactor(lam);
            mat->fForceFactor = lam;
        };
        p.setReduction = [&an](REAL F) { an.SetStrengthReductionFactor(F); };
        p.accept = [&an]() { QuietCout q; an.Accept(); };
        SetFreeEquations(an, p.free);
        const FactorResult f = ComputeFactor(p, res.neq, fFactor);
        res.value = f.value, res.steps = f.steps, res.end = f.end;
    }
    res.fmc = MaxYield(mat, d);
    res.ind = PlasticIndicator(cmesh, mat);
    if (!base.empty()) {
        WritePlasticPoints(cmesh, mat, base + "_pontos.csv");
        WritePointForces(cmesh, mat, base + "_forca.csv");
    }
    EraseMesh(cmesh);
    return res;
}

std::vector<LevelSummary> SlopeAnalysis::Run(FlatMesh m)
{
    for (int i = 0; i < fRef.nunif; i++) m = RefineUniform(m);
    const bool gim = fFactor == FactorMethod::GIMArcLength || fFactor == FactorMethod::GIMBisection;
    const std::string base = fPrefix.empty() ? "" : fPrefix + (fKind == MaterialKind::UP ? "_up" : "_unsat");
    std::ofstream resumo;
    if (!base.empty()) {
        resumo.open(base + ".csv");
        resumo << "nivel,triangulos,equacoes,livres,gamma,tempo_s,alvo\n";
    }
    const auto t0 = Clock::now();
    std::vector<LevelSummary> niveis;
    for (int k = 0;; k++) {
        const auto tk = Clock::now();
        LevelSummary ls;
        ls.level = k;
        ls.ntri = m.tri.size();
        for (int o : m.order) ls.orders[o]++;
        std::string ordens;
        for (const auto &[o, n] : ls.orders) ordens += " P" + std::to_string(o) + " " + std::to_string(n);
        const std::string header =
            "\n  nivel " + std::to_string(k) + " | " + std::to_string(m.tri.size()) + " triangulos (u:" + ordens + ")";
        const std::string nivel = base.empty() ? "" : base + "_nivel" + std::to_string(k);
        TPZGeoMesh *gmesh = ToGeoMesh(m);
        DuplicateBoundaryBC(gmesh);
        const LevelResult r = fKind == MaterialKind::UP ? SolveLevelUP(gmesh, m.order, header, nivel)
                                                        : SolveLevelUnsat(gmesh, m.order, header, nivel);
        delete gmesh;
        ls.neq = r.neq, ls.nfree = r.nfree, ls.value = r.value, ls.steps = r.steps, ls.seconds = Seconds(tk);
        niveis.push_back(ls);
        const REAL erro = fTarget > 0. ? 100. * (r.value / fTarget - 1.) : 0.;
        Print("    => %s = %.5f (erro %+.2f%%) | %d passos, fim: %s | f_MC max %.0e | %.1f s (total %.1f s)",
              gim ? "Gamma" : "FS", r.value, erro, r.steps, r.end.c_str(), r.fmc, Seconds(tk), Seconds(t0));
        if (!base.empty()) {
            std::vector<REAL> ordem(m.order.begin(), m.order.end());
            std::map<std::string, std::vector<REAL>> pontos;
            if (!r.pw.empty()) {
                pontos["p_w"] = r.pw;
                std::vector<REAL> u(m.nodes.size());
                for (size_t i = 0; i < u.size(); i++) u[i] = ExcessPorePressure(r.pw[i], m.nodes[i][1], fData);
                pontos["u_excesso"] = u;
            }
            WriteFlatVTK(m, {{"sqrtJ2_eps_p", std::vector<REAL>(r.ind.begin(), r.ind.begin() + m.tri.size())},
                             {"ordem_u", ordem}},
                         pontos, nivel + ".vtk");
            resumo << k << "," << m.tri.size() << "," << r.neq << "," << r.nfree << "," << r.value << ","
                   << Seconds(tk) << "," << fTarget << std::endl;
        }
        if (k >= fRef.nref) break;
        const REAL vmax = *std::max_element(r.ind.begin(), r.ind.begin() + m.tri.size());
        std::set<int64_t> marked;
        for (int64_t t = 0; t < (int64_t)m.tri.size(); t++)
            if (vmax > 0. && r.ind[t] >= fRef.mark * vmax) marked.insert(t);
        for (int64_t t : marked) m.order[t] = std::min(m.order[t] + 1, fRef.pmax); // p: os filhos herdam
        m = RefineMarked(m, marked);                                                // h
    }
    return niveis;
}

void PrintSummary(const std::string &title, const std::vector<SummaryRun> &runs, double total)
{
    int pmaxv = 2;
    size_t wnome = 8;
    for (const auto &r : runs) {
        wnome = std::max(wnome, r.name.size());
        for (const auto &l : r.levels) pmaxv = std::max(pmaxv, l.orders.rbegin()->first);
    }
    const std::string traco(wnome + 82 + 7 * (pmaxv - 1), '-'), duplo(traco.size(), '=');
    std::cout << "\n" << duplo << "\n RESUMO | " << title << "\n" << traco << "\n";
    char buf[160];
    std::snprintf(buf, sizeof buf, " %-*s %5s %8s %9s %8s", (int)wnome, "rodada", "nivel", "triang.", "equacoes",
                  "livres");
    std::cout << buf;
    for (int o = 2; o <= pmaxv; o++) std::cout << "     P" << o;
    std::cout << "      fator   alvo      erro  passos  tempo (s)\n";
    for (const auto &r : runs)
        for (const auto &l : r.levels) {
            std::snprintf(buf, sizeof buf, " %-*s %5d %8lld %9lld %8lld", (int)wnome, r.name.c_str(), l.level,
                          (long long)l.ntri, (long long)l.neq, (long long)l.nfree);
            std::cout << buf;
            for (int o = 2; o <= pmaxv; o++) {
                auto it = l.orders.find(o);
                std::snprintf(buf, sizeof buf, " %6d", it == l.orders.end() ? 0 : it->second);
                std::cout << buf;
            }
            Print(" %10.5f %6.3f %+8.2f%% %7d %10.1f", l.value, r.target, 100. * (l.value / r.target - 1.), l.steps,
                  l.seconds);
        }
    std::cout << traco << "\n";
    for (const auto &r : runs)
        Print(" final %-*s %.5f (alvo %.3f, %+.2f%%)", (int)wnome, r.name.c_str(), r.levels.back().value, r.target,
              100. * (r.levels.back().value / r.target - 1.));
    Print(" tempo total: %.1f s", total);
    std::cout << duplo << std::endl;
}

} // namespace slope
