// SlopeRandomField - talude de Rahardjo et al. (2007, Serie A, 45 graus) sob chuva com campos aleatorios de c', phi' e
// permeabilidade, gerados por Karhunen-Loeve com a estatistica de Vargas Ceron, Cecilio, Linn e Maghous (2025): campos
// lognormais e independentes, covariancia exponencial separavel exp(-|dx|/Lx - |dy|/Ly), k_h = alpha k_v. O fator de
// seguranca e calculado so por reducao de resistencia (SRM, com nova tentativa do F que falhou), nos instantes pedidos.
//
// Uma realizacao por execucao, na pasta corrente: REAL=<n> escolhe a realizacao (semente SEMENTE + n); CAMPOS=0 roda o
// caso deterministico (valores medios). O modelo (malha, solos, chuva, variantes) e o do driver de Rahardjo do projeto
// LearningMohrCoulomb2 (cluster/fonte/rahardjo.cpp), cujo main.cpp e incluido para reaproveitar malha, material e
// fatores de seguranca.
//
// Parametros (variaveis de ambiente) e saidas: ver README.md desta pasta.
#define main projeto_main_original
#include "../LearningMohrCoulomb2/cluster/fonte/main.cpp"
#undef main

#include "KLSeparavel.h"

#include <chrono>
#include <cstdlib>
#include <functional>
#include <sstream>

static double envd(const char *k, double d) { const char *v = std::getenv(k); return v ? std::atof(v) : d; }
static int envi(const char *k, int d) { const char *v = std::getenv(k); return v ? std::atoi(v) : d; }
static std::string envs(const char *k, const std::string &d) { const char *v = std::getenv(k); return v ? std::string(v) : d; }
static std::vector<double> envhoras(const char *k, const std::string &d) {
    std::vector<double> h;
    std::stringstream ss(envs(k, d));
    std::string t;
    while (std::getline(ss, t, ',')) if (!t.empty()) h.push_back(std::atof(t.c_str()));
    return h;
}

// FactorSRM do projeto + nova tentativa do F que falhou, agora com incremento < tol (copia de cluster/fonte/rahardjo.cpp)
static REAL FactorSRMRetry(TPZUnsatAnalysis &an, const TPZFMatrix<STATE> &x, REAL F0, REAL tol, const std::vector<REAL> &Fs)
{
    an.ResetMemory(); // mantem as propriedades de cada ponto (os campos)
    TPZFMatrix<STATE> &X = an.X();
    X = x;
    for (int64_t i = 0; i < an.NEq(); i++)
        if (!an.IsPressure(i)) X(i, 0) = 0.;
    an.Restore();
    const auto fixed = an.FrozenFixed(X);
    an.SetStrengthReductionFactor(F0);
    an.SetGravityScalesSuction(true);
    const int n = 4;
    for (int k = 1; k <= n; k++) {
        an.SetGravityFactor(REAL(k) / n);
        if (!an.Solve(0., fixed, {}, true).conv) return -1.; // sem equilibrio nem com c/F0, tan(phi)/F0: FS < F0
        an.Accept();
    }
    REAL Fok = F0, Ffail = -1.;
    size_t next = 0;
    while (true) {
        REAL F;
        bool retry = false;
        if (Ffail > 0.) {
            if (Ffail - Fok < tol) { F = Ffail; retry = true; }
            else F = 0.5 * (Fok + Ffail);
        } else {
            while (next < Fs.size() && Fs[next] <= Fok) next++;
            if (next == Fs.size()) break;
            F = Fs[next];
        }
        an.SetStrengthReductionFactor(F);
        const auto r = an.Solve(0., fixed, {}, true);
        if (r.conv) {
            an.Accept();
            Fok = F;
            if (retry) Ffail = -1.;
        } else {
            an.Restore();
            if (retry) break;
            Ffail = F;
        }
    }
    an.SetStrengthReductionFactor(1.);
    return Fok;
}

int main() {
    const auto inicio = std::chrono::steady_clock::now();
    const int porder = 2;
    const REAL hora = 3600.;
    // ---- chuva e tempo (Rahardjo et al. 2007, Serie A): 24 h de chuva; FS nos instantes FS_H
    const REAL I_mm_h = envd("I_MM_H", 36.), Ir = I_mm_h / (1000. * hora);
    const REAL Tchuva = envd("TCHUVA_H", 24.) * hora, Dt = envd("DT_H", 0.1) * hora;
    const std::vector<double> fsH = envhoras("FS_H", "0,12,24,36,48");
    const REAL Tf = std::max(envd("TF_H", 0.), *std::max_element(fsH.begin(), fsH.end())) * hora;
    const REAL xl = envd("XL", 0.), xr = envd("XR", 70.);
    const REAL Hw = envd("HW", 5.), beta = envd("BETA_DEG", 7.) * M_PI / 180., sucMax = envd("SUCMAX", 75.);
    const REAL xtoe = envd("XTOE", 40.), ytoe = envd("YTOE", 20.);
    const REAL hl = envd("HL", ytoe - Hw + (xtoe - xl) * std::tan(beta)), hr = envd("HR", ytoe - Hw + (xtoe - xr) * std::tan(beta));
    const int verticalRain = envi("CHUVA_VERTICAL", 0);
    const REAL fsTol = envd("FS_TOL", 0.01), F0 = envd("FS_F0", 0.6);
    const int maxCortes = envi("MAX_CORTES", 4); // passo que nao converge: ate 2^MAX_CORTES subpassos

    // ---- solo (Rahardjo et al. 2007, Tabela 1): c' = 10 kPa, phi' = phi^b = 26 graus, gamma = 20 kN/m3, n = 0.45;
    // Fredlund-Xing com a = 10/50/100 kPa, m = n = 1; k_w = k_s Theta^4 com k_s = 1e-4/1e-5/1e-6 m/s
    const std::string solo = envs("SOLO", "f50");
    const REAL aFX = solo == "f10" ? 10. : solo == "f50" ? 50. : solo == "f100" ? 100. : -1.;
    if (aFX < 0.) { std::cerr << "SOLO deve ser f10, f50 ou f100\n"; return 2; }
    const std::string variante = envs("VARIANTE", "iii"); // i: completo; ii: ar passivo; iii: esqueleto rigido (SEEP/W)
    UnsatParams p;
    p.E = 2.e4, p.nu = 0.3, p.c = 10., p.phi = 26. * M_PI / 180., p.gamma = 20.;
    p.n = 0.45, p.ks = solo == "f10" ? 1.e-4 : solo == "f50" ? 1.e-5 : 1.e-6;
    p.hydro.fRetention = TPZUnsaturatedHydraulics::EFredlundXing;
    p.hydro.fKrw = TPZUnsaturatedHydraulics::EPowerLaw;
    p.hydro.fKra = TPZUnsaturatedHydraulics::EBrooksCorey;
    p.hydro.fA = aFX, p.hydro.fN = 1., p.hydro.fM = 1., p.hydro.fPower = 4.;
    p.bishop = TUnsat::EChiConstant, p.chi = 1.; // phi^b = phi' em cada ponto
    p.coupled = variante == "i" || variante == "ii";
    p.E = envd("E", p.E); p.nu = envd("NU", p.nu); p.c = envd("COES", p.c); p.phi = envd("PHI_DEG", p.phi * 180. / M_PI) * M_PI / 180.;
    p.gamma = envd("GAMMA", p.gamma); p.ks = envd("KS", p.ks); p.n = envd("POROS", p.n);

    // ---- estatistica dos campos (Vargas Ceron et al. 2025, Tabela 2): CoV(c') 30%, CoV(phi') 10%, CoV(k_v) 60%,
    // Lx = 20 m, Ly = 2 m, alpha = k_h/k_v = 1; medias de c', phi' e k_s do solo de Rahardjo; M = 2000 termos de KL
    const int campos = envi("CAMPOS", 1), real = envi("REAL", 0);
    const uint64_t semente = (uint64_t)envd("SEMENTE", 20250416.);
    const REAL covC = envd("COV_C", 0.3), covPhi = envd("COV_PHI", 0.1), covK = envd("COV_K", 0.6);
    const REAL Lx = envd("LX", 20.), Ly = envd("LY", 2.), alfa = envd("ALFA", 1.);
    const int Mkl = envi("KL_M", 2000);
    const REAL fracKL = envd("KL_FRAC", 1.);

    TPZGeoMesh *gmesh = ReadGmshMesh(envs("MESH", "malhas/rahardjo_a45.msh"), {{"Solo", matSolo}},
        {{"Base", bcBase}, {"Direita", bcDireita}, {"Superficie", bcSuperficie}, {"Esquerda", bcEsquerda}});
    TPZMultiphysicsCompMesh *cmesh = CreateSlopeCompMesh(gmesh, porder, p);
    TPZUnsatAnalysis::SetUseSkyline(envs("SOLVER", "pardiso") == "skyline");
    auto *an = new TPZUnsatAnalysis(cmesh, {{bcBase, {true, true}}, {bcDireita, {true, false}}, {bcEsquerda, {true, false}}});
    an->SetTolerances(envd("TOLU", 1.e-7), envd("TOLP", 1.e-6), envi("MAXIT", 40), envi("MAXIT_FROZEN", 25));
    auto *solomat = dynamic_cast<TUnsat *>(cmesh->FindMaterial(matSolo));
    if (!solomat) DebugStop();
    // tensor de condutividade medio: k_h = alpha k_v na horizontal, k_v = k_s na vertical
    solomat->SetAnisotropicConductivity(alfa * p.ks, p.ks, 0.);
    const auto borda = an->BoundaryNodes(bcSuperficie, verticalRain != 0);
    if (borda.empty()) DebugStop();

    // ---- pontos de integracao do solo: indice de memoria e coordenadas
    std::vector<int64_t> ips;
    std::vector<double> xs, ys;
    for (TPZCompEl *cel : cmesh->ElementVec()) {
        if (!cel || !cel->Reference() || cel->Reference()->MaterialId() != matSolo) continue;
        TPZManVector<int64_t, 10> idx;
        cel->GetMemoryIndices(idx);
        const TPZIntPoints &rule = cel->GetIntegrationRule();
        if (idx.size() != rule.NPoints()) DebugStop();
        TPZManVector<REAL, 3> qsi(2, 0.), x(3, 0.);
        REAL w;
        for (int i = 0; i < rule.NPoints(); ++i) {
            rule.Point(i, qsi, w);
            cel->Reference()->X(qsi, x);
            ips.push_back(idx[i]);
            xs.push_back(x[0]);
            ys.push_back(x[1]);
        }
    }

    // ---- campos de Karhunen-Loeve no retangulo que contem a malha, aplicados a cada ponto
    std::ofstream info("realizacao.txt");
    info << std::setprecision(10) << "solo=" << solo << " I_mm_h=" << I_mm_h << " variante=" << variante
         << " solver=" << (TPZUnsatAnalysis::UsesPardiso() ? "pardiso" : "skyline") << " neq=" << an->NEq()
         << " pontos=" << ips.size() << " real=" << real << " campos=" << campos;
    if (campos) {
        REAL x0 = 1.e30, x1 = -1.e30, y0 = 1.e30, y1 = -1.e30;
        for (int64_t i = 0; i < gmesh->NNodes(); ++i) {
            TPZManVector<REAL, 3> c(3);
            gmesh->NodeVec()[i].GetCoordinates(c);
            x0 = std::min(x0, c[0]); x1 = std::max(x1, c[0]); y0 = std::min(y0, c[1]); y1 = std::max(y1, c[1]);
        }
        const klsep::KL2D kl(x0, y0, x1 - x0, y1 - y0, Lx, Ly, fracKL, Mkl);
        const auto P = kl.PreparaPontos(xs, ys);
        // tres sequencias independentes por realizacao (c', phi', k): sementes derivadas de (SEMENTE, REAL, campo)
        auto xi = [&](int campo) { return klsep::Normal(semente * 0x9E3779B97F4A7C15ULL + 3ULL * real + campo).Vetor(kl.M()); };
        const auto zc = kl.Campo(P, xi(0)), zphi = kl.Campo(P, xi(1)), zk = kl.Campo(P, xi(2));
        const klsep::Lognormal lc(p.c, covC), lphi(p.phi * 180. / M_PI, covPhi), lk(1., covK);
        auto mem = solomat->GetMemory();
        std::ofstream fc;
        if (envi("CAMPOS_CSV", 0)) { fc.open("campos.csv"); fc << "ip,x,y,c,phi_deg,k_mult\n" << std::setprecision(8); }
        double sc = 0., sphi = 0., sk = 0.;
        for (size_t i = 0; i < ips.size(); ++i) {
            const REAL c = lc(zc[i]), phi = lphi(zphi[i]) * M_PI / 180., km = lk(zk[i]);
            TPZUnsaturatedMem &m = (*mem)[ips[i]];
            m.m_elastoplastic_state.fmatprop = {c, phi, phi}; // {c, phi, psi}, psi = phi (associado)
            m.fKMult = km;
            sc += c; sphi += phi * 180. / M_PI; sk += km;
            if (fc.is_open()) fc << ips[i] << ',' << xs[i] << ',' << ys[i] << ',' << c << ',' << phi * 180. / M_PI << ',' << km << '\n';
        }
        const double np = double(ips.size());
        info << " KL_M=" << kl.M() << " KL_variancia=" << kl.fracao << " Lx=" << Lx << " Ly=" << Ly << " alfa=" << alfa
             << " cov_c=" << covC << " cov_phi=" << covPhi << " cov_k=" << covK << " media_c=" << sc / np
             << " media_phi=" << sphi / np << " media_kmult=" << sk / np;
    }
    info << '\n';
    info.close();

    // ---- estado inicial (como no driver de Rahardjo): hidrostatico com succao limitada; equilibrio com a resistencia real
    TPZFMatrix<STATE> x(cmesh->NEquations(), 1, 0.);
    for (const auto &eq : an->WaterEq()) {
        const REAL xx = an->Coord(eq.first, 0), y = an->Coord(eq.first, 1);
        const REAL hf = hl + (hr - hl) * (xx - xl) / (xr - xl);
        x(eq.second, 0) = std::max(p.gw * (hf - y), -sucMax);
    }
    an->X() = x;
    an->Restore();
    {
        const auto inicial = an->FrozenFixed(x);
        an->SetStrengthReductionFactor(1.);
        an->SetGravityScalesSuction(true);
        bool ok = true;
        for (int k = 1; k <= 4 && ok; k++) {
            an->SetGravityFactor(REAL(k) / 4.);
            ok = an->Solve(0., inicial, {}, true).conv;
            if (ok) an->Accept();
        }
        if (ok) an->ResetDisplacements();
        else std::cerr << "estado inicial sem equilibrio com a resistencia desta realizacao\n";
    }

    // ---- contornos da percolacao (como no driver de Rahardjo)
    auto contorno = variante == "iii" ? an->AllDisplacements(an->X()) : an->DisplacementConstraints();
    for (const int id : {bcEsquerda, bcDireita}) {
        const REAL h = id == bcEsquerda ? hl : hr;
        for (const auto &b : an->BoundaryNodes(id)) {
            const REAL y = an->Coord(b.first, 1);
            if (y <= h) contorno[an->WaterEq().at(b.first)] = p.gw * (h - y);
        }
    }
    for (const auto &b : borda) contorno[an->AirEq().at(b.first)] = 0.;
    if (variante != "i") for (const auto &e : an->AirEq()) contorno[e.second] = 0.;

    // ---- um passo de chuva de t0 a t1, com a face de percolacao (o que estiver com p_w > 0 na superficie fica em p_w = 0
    // enquanto a infiltracao exigida for menor que a chuva); em caso de falha, volta ao estado de t0
    std::set<int64_t> saida;
    const REAL tolP = 1.e-6, tolQ = 1.e-10;
    int64_t nPassos = 0, nCortes = 0, newtonTotal = 0;
    auto passo = [&](REAL t0, REAL t1) -> bool {
        const REAL dt = t1 - t0, chuva = Ir * std::max(REAL(0), std::min(t1, Tchuva) - t0) / dt;
        auto tentativa = saida;
        for (const auto &b : borda)
            if (an->X().GetVal(an->WaterEq().at(b.first), 0) > tolP) tentativa.insert(b.first);
        for (int it = 0; it < 40; it++) {
            auto fixed = contorno;
            std::map<int64_t, REAL> fq;
            for (const auto &b : borda) {
                const int64_t eq = an->WaterEq().at(b.first);
                if (tentativa.count(b.first)) fixed[eq] = 0.;
                else fq[eq] = chuva * b.second;
            }
            const auto sr = an->Solve(dt, fixed, fq, false);
            newtonTotal += sr.iter;
            if (!sr.conv) break;
            const auto R = an->FullResidual();
            const TPZFMatrix<STATE> sol = cmesh->Solution();
            auto nova = tentativa;
            for (const auto &b : borda) {
                const int64_t eq = an->WaterEq().at(b.first);
                if (tentativa.count(b.first)) {
                    if (R.GetVal(eq, 0) / (dt * b.second) > chuva + tolQ) nova.erase(b.first);
                } else if (sol.GetVal(eq, 0) > tolP) nova.insert(b.first);
            }
            if (nova == tentativa) {
                an->Accept();
                saida = tentativa;
                nPassos++;
                return true;
            }
            tentativa = nova;
        }
        an->Restore();
        return false;
    };
    // avanca de t0 a t1, dividindo ao meio o passo que nao converge (ate maxCortes niveis)
    std::function<bool(REAL, REAL, int)> avanca = [&](REAL t0, REAL t1, int nivel) -> bool {
        if (passo(t0, t1)) return true;
        if (nivel >= maxCortes) return false;
        nCortes++;
        const REAL tm = 0.5 * (t0 + t1);
        return avanca(t0, tm, nivel + 1) && avanca(tm, t1, nivel + 1);
    };

    // ---- transiente com FS por SRM nos instantes pedidos
    std::vector<REAL> Fdenso;
    for (REAL f = 0.7; f < 3.05; f += 0.1) Fdenso.push_back(std::round(f * 10.) / 10.);
    Fdenso.push_back(3.5);
    Fdenso.push_back(4.);
    std::ofstream fs("fs.csv");
    fs << "t_h,FS_SRM,seg\n" << std::setprecision(10);
    const int nsteps = int(std::lround(Tf / Dt));
    std::set<int> fsSteps;
    for (double h : fsH) fsSteps.insert(int(std::lround(h * hora / Dt)));
    REAL t = 0., fsMin = 1.e30, tMin = -1.;
    std::string situacao = "ok";
    for (int step = 0; step <= nsteps; step++) {
        if (step > 0) {
            const REAL tn = std::min(step * Dt, Tf);
            if (!avanca(t, tn, 0)) {
                situacao = "transiente_nao_convergiu_t=" + std::to_string(t / hora);
                break;
            }
            t = tn;
        }
        if (fsSteps.count(step)) {
            const auto t0 = std::chrono::steady_clock::now();
            const auto estado = an->X();
            const auto memoria = an->SaveMemory();
            const REAL dtAnterior = solomat->TimeStep();
            const REAL F = FactorSRMRetry(*an, estado, F0, fsTol, Fdenso);
            an->SetStrengthReductionFactor(1.);
            an->SetGravityFactor(1.);
            an->SetGravityScalesSuction(true);
            an->RestoreMemory(memoria);
            an->X() = estado;
            an->Restore();
            solomat->SetTimeStep(dtAnterior);
            solomat->SetFrozenPressures(false);
            const double seg = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
            fs << t / hora << ',' << F << ',' << seg << '\n';
            fs.flush();
            std::cout << "t = " << t / hora << " h: FS = " << F << " (" << seg << " s)" << std::endl;
            if (F < fsMin) { fsMin = F; tMin = t / hora; }
        }
    }
    const double total = std::chrono::duration<double>(std::chrono::steady_clock::now() - inicio).count();
    std::ofstream res("resultado.txt");
    res << std::setprecision(10) << "real=" << real << " FS_min=" << fsMin << " t_min_h=" << tMin << " situacao=" << situacao
        << " passos=" << nPassos << " cortes=" << nCortes << " newton=" << newtonTotal << " seg=" << total << '\n';
    std::cout << "real " << real << ": FS_min = " << fsMin << " em t = " << tMin << " h, " << situacao << ", " << total << " s"
              << std::endl;
    delete an;
    DeleteCompMesh(cmesh);
    delete gmesh;
    return situacao == "ok" ? 0 : 1;
}
