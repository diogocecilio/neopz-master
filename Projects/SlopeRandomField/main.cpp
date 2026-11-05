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
// Saidas em CSV, para o pos-processamento: resumo.csv (uma linha com tudo da realizacao), fs.csv (cada instante do FS,
// com o mecanismo e o estado hidraulico), fs_tentativas.csv (cada solucao da busca do SRM), passos.csv (cada passo do
// transiente, com o balanco de agua), pw.csv (p_w nos nos nos instantes do FS), mecanismo.csv (pontos que plastificaram
// no ultimo estado convergido do SRM) e campos.csv (c', phi' e k em cada ponto de integracao). Uma execucao interrompida
// retoma do checkpoint.bin. Parametros e colunas: README.md desta pasta.
#define main projeto_main_original
#include "../LearningMohrCoulomb2/cluster/fonte/main.cpp"
#undef main

#include "KLSeparavel.h"
#include "TPZBFileStream.h"

#include <array>
#include <chrono>
#include <cstdlib>
#include <filesystem>
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

using Relogio = std::chrono::steady_clock;
static double Segundos(Relogio::time_point t0) { return std::chrono::duration<double>(Relogio::now() - t0).count(); }

// pico da memoria residente do processo (VmHWM), em MB; -1 fora do Linux
static double PicoMemoriaMB() {
    std::ifstream f("/proc/self/status");
    std::string l;
    while (std::getline(f, l))
        if (l.rfind("VmHWM:", 0) == 0) return std::atof(l.c_str() + 6) / 1024.;
    return -1.;
}

static std::string Num(double v) { std::ostringstream o; o << std::setprecision(10) << v; return o.str(); }

static std::vector<std::string> SeparaCsv(const std::string &l) {
    std::vector<std::string> v;
    std::stringstream ss(l);
    std::string x;
    while (std::getline(ss, x, ',')) v.push_back(x);
    if (!l.empty() && l.back() == ',') v.push_back("");
    return v;
}

// uma solucao da busca do SRM: um dos 4 passos de carga da gravidade com c'/F0, tan(phi')/F0 ("gravidade") ou um F
// tentado depois ("srm")
struct TentativaSRM {
    const char *fase;
    REAL g, F;
    bool conv, retry;
    int iter;
    double seg;
};

// FactorSRM do projeto + nova tentativa do F que falhou, agora com incremento < tol (como em cluster/fonte/rahardjo.cpp).
// Fparada: o F em que o FactorSRM do projeto pararia, sem a nova tentativa; tent: cada solucao tentada
static REAL FactorSRMRetry(TPZUnsatAnalysis &an, const TPZFMatrix<STATE> &x, REAL F0, REAL tol, const std::vector<REAL> &Fs,
                           REAL &Fparada, std::vector<TentativaSRM> &tent)
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
    Fparada = -1.;
    const int n = 4;
    for (int k = 1; k <= n; k++) {
        an.SetGravityFactor(REAL(k) / n);
        const auto t0 = Relogio::now();
        const auto r = an.Solve(0., fixed, {}, true);
        tent.push_back({"gravidade", REAL(k) / n, F0, r.conv, false, r.iter, Segundos(t0)});
        if (!r.conv) return -1.; // sem equilibrio nem com c'/F0, tan(phi')/F0: FS < F0
        an.Accept();
    }
    REAL Fok = F0, Ffail = -1.;
    size_t next = 0;
    while (true) {
        REAL F;
        bool retry = false;
        if (Ffail > 0.) {
            if (Ffail - Fok < tol) {
                if (Fparada < 0.) Fparada = Fok;
                F = Ffail;
                retry = true;
            } else F = 0.5 * (Fok + Ffail);
        } else {
            while (next < Fs.size() && Fs[next] <= Fok) next++;
            if (next == Fs.size()) break;
            F = Fs[next];
        }
        an.SetStrengthReductionFactor(F);
        const auto t0 = Relogio::now();
        const auto r = an.Solve(0., fixed, {}, true);
        tent.push_back({"srm", 1., F, r.conv, retry, r.iter, Segundos(t0)});
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
    if (Fparada < 0.) Fparada = Fok;
    an.SetStrengthReductionFactor(1.);
    return Fok;
}

int main() {
    const auto inicio = Relogio::now();
    const int porder = 2;
    const REAL hora = 3600.;
    // ---- chuva e tempo (Rahardjo et al. 2007, Serie A): 24 h de chuva; FS nos instantes FS_H
    const REAL I_mm_h = envd("I_MM_H", 36.), Ir = I_mm_h / (1000. * hora);
    const REAL Tchuva = envd("TCHUVA_H", 24.) * hora, Dt = envd("DT_H", 0.1) * hora;
    const std::vector<double> fsH = envhoras("FS_H", "0,12,24,36,48");
    if (fsH.empty()) { std::cerr << "FS_H sem instantes\n"; return 2; }
    const REAL Tf = std::max(envd("TF_H", 0.), *std::max_element(fsH.begin(), fsH.end())) * hora;
    const REAL xl = envd("XL", 0.), xr = envd("XR", 70.);
    const REAL Hw = envd("HW", 5.), beta = envd("BETA_DEG", 7.) * M_PI / 180., sucMax = envd("SUCMAX", 75.);
    const REAL xtoe = envd("XTOE", 40.), ytoe = envd("YTOE", 20.);
    const REAL hl = envd("HL", ytoe - Hw + (xtoe - xl) * std::tan(beta)), hr = envd("HR", ytoe - Hw + (xtoe - xr) * std::tan(beta));
    const int verticalRain = envi("CHUVA_VERTICAL", 0);
    const REAL fsTol = envd("FS_TOL", 0.01), F0 = envd("FS_F0", 0.6);
    const int maxCortes = envi("MAX_CORTES", 4);      // passo que nao converge: ate 2^MAX_CORTES subpassos
    const double ckMin = envd("CHECKPOINT_MIN", 5.); // checkpoint a cada CHECKPOINT_MIN minutos e depois de cada FS
    // CSV opcionais (1 = grava): p_w nos nos, pontos do mecanismo e campos nos pontos de integracao
    const int pwCsv = envi("PW_CSV", 1), mecCsv = envi("MEC_CSV", 1), camposCsv = envi("CAMPOS_CSV", 1);
    // mecanismo.csv: pontos com sqrt(J2(eps_p)) >= MEC_LIMIAR * maximo; banda de ruptura (fs.csv): >= BANDA_LIMIAR * maximo
    const REAL mecLimiar = envd("MEC_LIMIAR", 0.01), bandaLimiar = envd("BANDA_LIMIAR", 0.1);

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

    // ---- pontos de integracao do solo: indice de memoria, elemento, coordenadas e area (peso vezes det J)
    std::vector<int64_t> ips, ipel;
    std::vector<double> xs, ys, area;
    for (TPZCompEl *cel : cmesh->ElementVec()) {
        if (!cel || !cel->Reference() || cel->Reference()->MaterialId() != matSolo) continue;
        TPZManVector<int64_t, 10> idx;
        cel->GetMemoryIndices(idx);
        const TPZIntPoints &rule = cel->GetIntegrationRule();
        if (idx.size() != rule.NPoints()) DebugStop();
        TPZManVector<REAL, 3> qsi(2, 0.), x(3, 0.);
        TPZFNMatrix<9, REAL> jac, axes, jacinv;
        REAL w, detjac;
        for (int i = 0; i < rule.NPoints(); ++i) {
            rule.Point(i, qsi, w);
            cel->Reference()->X(qsi, x);
            cel->Reference()->Jacobian(qsi, jac, axes, detjac, jacinv);
            ips.push_back(idx[i]);
            ipel.push_back(cel->Reference()->Index());
            xs.push_back(x[0]);
            ys.push_back(x[1]);
            area.push_back(w * std::fabs(detjac));
        }
    }
    auto mem = solomat->GetMemory();
    // propriedades do ponto i (os valores do material se o ponto nao tem campo: CAMPOS=0); phi em graus
    auto cPonto = [&](size_t i) {
        const auto &mp = (*mem)[ips[i]].m_elastoplastic_state.fmatprop;
        return mp.size() >= 2 ? mp[0] : p.c;
    };
    auto phiPonto = [&](size_t i) {
        const auto &mp = (*mem)[ips[i]].m_elastoplastic_state.fmatprop;
        return (mp.size() >= 2 ? mp[1] : p.phi) * 180. / M_PI;
    };
    auto kPonto = [&](size_t i) { return (*mem)[ips[i]].fKMult; };

    // ---- campos de Karhunen-Loeve no retangulo que contem a malha, aplicados a cada ponto
    std::ofstream info("realizacao.txt");
    info << std::setprecision(10) << "solo=" << solo << " I_mm_h=" << I_mm_h << " variante=" << variante
         << " solver=" << (TPZUnsatAnalysis::UsesPardiso() ? "pardiso" : "skyline") << " neq=" << an->NEq()
         << " pontos=" << ips.size() << " real=" << real << " campos=" << campos;
    int Musado = 0;
    double fracUsada = 0.;
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
        for (size_t i = 0; i < ips.size(); ++i) {
            const REAL c = lc(zc[i]), phi = lphi(zphi[i]) * M_PI / 180., km = lk(zk[i]);
            TPZUnsaturatedMem &m = (*mem)[ips[i]];
            m.m_elastoplastic_state.fmatprop = {c, phi, phi}; // {c, phi, psi}, psi = phi (associado)
            m.fKMult = km;
        }
        Musado = kl.M();
        fracUsada = kl.fracao;
        if (camposCsv) {
            std::ofstream fc("campos.csv");
            fc << "ip,x,y,area,c,phi_deg,k_mult\n" << std::setprecision(7);
            for (size_t i = 0; i < ips.size(); ++i)
                fc << ips[i] << ',' << xs[i] << ',' << ys[i] << ',' << area[i] << ',' << cPonto(i) << ',' << phiPonto(i) << ','
                   << kPonto(i) << '\n';
        }
        if (envi("CAMPOS_VTK", 0)) {
            // um valor por triangulo, a media dos seus pontos de integracao, sobre a malha geometrica (cantos)
            std::map<int64_t, std::array<double, 4>> soma; // c, phi, k_mult, numero de pontos
            for (size_t i = 0; i < ips.size(); ++i) {
                auto &a = soma[ipel[i]];
                a[0] += cPonto(i); a[1] += phiPonto(i); a[2] += kPonto(i); a[3] += 1.;
            }
            std::ofstream v("campos.vtk");
            v << "# vtk DataFile Version 3.0\nSlopeRandomField realizacao " << real << "\nASCII\nDATASET UNSTRUCTURED_GRID\n";
            v << std::setprecision(8) << "POINTS " << gmesh->NNodes() << " double\n";
            for (int64_t i = 0; i < gmesh->NNodes(); ++i) {
                TPZManVector<REAL, 3> c(3, 0.);
                gmesh->NodeVec()[i].GetCoordinates(c);
                v << c[0] << ' ' << c[1] << " 0\n";
            }
            v << "CELLS " << soma.size() << ' ' << 4 * soma.size() << '\n';
            for (const auto &e : soma) {
                TPZGeoEl *g = gmesh->Element(e.first);
                v << "3 " << g->NodeIndex(0) << ' ' << g->NodeIndex(1) << ' ' << g->NodeIndex(2) << '\n';
            }
            v << "CELL_TYPES " << soma.size() << '\n';
            for (size_t i = 0; i < soma.size(); ++i) v << "5\n";
            v << "CELL_DATA " << soma.size() << '\n';
            auto escalar = [&](const char *nome, auto f) {
                v << "SCALARS " << nome << " double 1\nLOOKUP_TABLE default\n";
                for (const auto &e : soma) v << f(e.second) << '\n';
            };
            escalar("coesao_kPa", [](const std::array<double, 4> &a) { return a[0] / a[3]; });
            escalar("atrito_graus", [](const std::array<double, 4> &a) { return a[1] / a[3]; });
            escalar("k_v_m_s", [&](const std::array<double, 4> &a) { return p.ks * a[2] / a[3]; });
            escalar("log10_k_v", [&](const std::array<double, 4> &a) { return std::log10(p.ks * a[2] / a[3]); });
            escalar("k_h_m_s", [&](const std::array<double, 4> &a) { return alfa * p.ks * a[2] / a[3]; });
        }
    }
    // estatistica de cada campo no solo, ponderada pela area: media, desvio, minimo e maximo
    struct Estat { double media = 0., dp = 0., min = 1.e300, max = -1.e300; };
    auto estatistica = [&](auto f) {
        Estat e;
        double sa = 0., s1 = 0., s2 = 0.;
        for (size_t i = 0; i < ips.size(); ++i) {
            const double v = f(i), a = area[i];
            sa += a; s1 += a * v; s2 += a * v * v;
            e.min = std::min(e.min, v); e.max = std::max(e.max, v);
        }
        e.media = s1 / sa;
        e.dp = std::sqrt(std::max(0., s2 / sa - e.media * e.media));
        return e;
    };
    const Estat ec = estatistica(cPonto), ephi = estatistica(phiPonto), ek = estatistica(kPonto);
    double areaSolo = 0.;
    for (double a : area) areaSolo += a;
    if (campos)
        info << " KL_M=" << Musado << " KL_variancia=" << fracUsada << " Lx=" << Lx << " Ly=" << Ly << " alfa=" << alfa
             << " cov_c=" << covC << " cov_phi=" << covPhi << " cov_k=" << covK;
    info << " media_c=" << ec.media << " media_phi=" << ephi.media << " media_kmult=" << ek.media << '\n';
    info.close();
    if (envi("SO_CAMPOS", 0)) { // so os campos (campos.csv/campos.vtk), sem simulacao: a mesma SEMENTE e REAL refazem o campo
        delete an;
        DeleteCompMesh(cmesh);
        delete gmesh;
        return 0;
    }

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

    // ---- contornos da percolacao (como no driver de Rahardjo); eqEsq/eqDir: nos com p_w imposto nas laterais
    auto contorno = variante == "iii" ? an->AllDisplacements(an->X()) : an->DisplacementConstraints();
    std::set<int64_t> eqEsq, eqDir;
    for (const int id : {bcEsquerda, bcDireita}) {
        const REAL h = id == bcEsquerda ? hl : hr;
        for (const auto &b : an->BoundaryNodes(id)) {
            const REAL y = an->Coord(b.first, 1);
            if (y <= h) {
                contorno[an->WaterEq().at(b.first)] = p.gw * (h - y);
                (id == bcEsquerda ? eqEsq : eqDir).insert(an->WaterEq().at(b.first));
            }
        }
    }
    for (const auto &b : borda) contorno[an->AirEq().at(b.first)] = 0.;
    if (variante != "i") for (const auto &e : an->AirEq()) contorno[e.second] = 0.;

    // volume de agua no solo (m3 por metro de talude), com a massa concentrada de cada vertice (area/3 dos triangulos);
    // triangulos.csv: os nos de cada triangulo (os numeros de no de pw.csv), para desenhar p_w sobre a malha
    std::map<int64_t, REAL> massa;
    std::ofstream ftri("triangulos.csv");
    ftri << "n0,n1,n2\n";
    for (int64_t i = 0; i < gmesh->NElements(); i++) {
        TPZGeoEl *g = gmesh->Element(i);
        if (!g || g->Dimension() != 2 || g->MaterialId() != matSolo) continue;
        ftri << g->NodeIndex(0) << ',' << g->NodeIndex(1) << ',' << g->NodeIndex(2) << '\n';
        const REAL x0 = an->Coord(g->NodeIndex(0), 0), y0 = an->Coord(g->NodeIndex(0), 1);
        const REAL x1 = an->Coord(g->NodeIndex(1), 0), y1 = an->Coord(g->NodeIndex(1), 1);
        const REAL x2 = an->Coord(g->NodeIndex(2), 0), y2 = an->Coord(g->NodeIndex(2), 1);
        const REAL A = 0.5 * std::fabs((x1 - x0) * (y2 - y0) - (x2 - x0) * (y1 - y0));
        for (int k = 0; k < 3; k++) massa[g->NodeIndex(k)] += A / 3.;
    }
    auto VolumeAgua = [&](const TPZFMatrix<STATE> &xv) {
        REAL v = 0.;
        for (const auto &m : massa) {
            const REAL pw = xv.GetVal(an->WaterEq().at(m.first), 0), pa = xv.GetVal(an->AirEq().at(m.first), 0);
            v += m.second * p.n * solomat->Hydraulics().Evaluate(pw, pa).S;
        }
        return v;
    };

    // ---- arquivos CSV: abertos depois de ler o checkpoint (na retomada, cortados no tamanho que tinham nele)
    std::ofstream fpas, ffs, ften, fpw, fmec;
    struct Csv { const char *nome; std::ofstream *o; const char *cab; bool ativo; int prec; };
    const std::vector<Csv> csvs = {
        {"passos.csv", &fpas,
         "passo,t_h,dt_h,nivel,chuva_mm_h,chuva_m3,inf_fluxo_m3,inf_face_m3,esq_m3,dir_m3,dV_m3,residuo_m3,inf_acum_m3,"
         "V_agua_m3,face_nos,its_face,newton,normp,eu,pw_min_kPa,pw_max_kPa,seg", true, 10},
        {"fs.csv", &ffs,
         "t_h,FS_SRM,FS_parada,tentativas,newton_fs,u_max_m,pl_pts,epsp_max,banda_pts,banda_area_m2,banda_xc,banda_yc,"
         "banda_xmin,banda_xmax,banda_ymin,banda_ymax,banda_c,banda_phi_deg,banda_k_mult,inf_acum_m3,V_agua_m3,pw_max_kPa,"
         "face_nos,seg", true, 10},
        {"fs_tentativas.csv", &ften, "t_h,ordem,fase,g,F,convergiu,newton,nova_tentativa,seg", true, 10},
        {"pw.csv", &fpw, "t_h,no,x,y,pw_kPa,pa_kPa,S,superficie,face", pwCsv != 0, 7},
        {"mecanismo.csv", &fmec, "t_h,FS,ip,x,y,epsp_raizJ2,epsp_xx,epsp_yy,epsp_xy,plastico,c,phi_deg,k_mult", mecCsv != 0, 7}};

    // ---- um passo de chuva de t0 a t1, com a face de percolacao (o que estiver com p_w > 0 na superficie fica em p_w = 0
    // enquanto a infiltracao exigida for menor que a chuva); em caso de falha, volta ao estado de t0. Cada passo aceito
    // vai para passos.csv com o balanco de agua (volumes por metro de talude no passo): R_w e a entrada de agua nos nos
    // com p_w imposto e a soma de R_w em todos os nos e a variacao do armazenamento (dV)
    std::set<int64_t> saida;
    const REAL tolP = 1.e-6, tolQ = 1.e-10;
    int64_t nPassos = 0, nCortes = 0, newtonTotal = 0, retomadas = 0;
    double infAcum = 0.; // infiltracao liquida acumulada pela superficie (m3 por metro)
    int stepAtual = 0;   // passo de Dt a que o subpasso pertence
    auto passo = [&](REAL t0, REAL t1, int nivel) -> bool {
        const auto relogio0 = Relogio::now();
        const REAL dt = t1 - t0, chuva = Ir * std::max(REAL(0), std::min(t1, Tchuva) - t0) / dt;
        auto tentativa = saida;
        for (const auto &b : borda)
            if (an->X().GetVal(an->WaterEq().at(b.first), 0) > tolP) tentativa.insert(b.first);
        int newtonPasso = 0;
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
            newtonPasso += sr.iter;
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
                REAL somaR = 0., chuvaPot = 0., infFluxo = 0., infFace = 0., qEsq = 0., qDir = 0.;
                REAL pwmin = 1.e300, pwmax = -1.e300;
                for (const auto &e : an->WaterEq()) {
                    somaR += R.GetVal(e.second, 0);
                    const REAL pw = sol.GetVal(e.second, 0);
                    pwmin = std::min(pwmin, pw);
                    pwmax = std::max(pwmax, pw);
                }
                for (const auto &b : borda) {
                    chuvaPot += chuva * b.second * dt;
                    if (tentativa.count(b.first)) infFace += R.GetVal(an->WaterEq().at(b.first), 0);
                    else infFluxo += chuva * b.second * dt;
                }
                for (auto e : eqEsq) qEsq += R.GetVal(e, 0);
                for (auto e : eqDir) qDir += R.GetVal(e, 0);
                const REAL V = VolumeAgua(sol);
                an->Accept();
                saida = tentativa;
                nPassos++;
                infAcum += infFluxo + infFace;
                fpas << stepAtual << ',' << t1 / hora << ',' << dt / hora << ',' << nivel << ',' << chuva * 1000. * hora << ','
                     << chuvaPot << ',' << infFluxo << ',' << infFace << ',' << qEsq << ',' << qDir << ',' << somaR << ','
                     << somaR - infFluxo - infFace - qEsq - qDir << ',' << infAcum << ',' << V << ',' << saida.size() << ','
                     << it + 1 << ',' << newtonPasso << ',' << sr.normp << ',' << sr.eu << ',' << pwmin << ',' << pwmax << ','
                     << Segundos(relogio0) << '\n';
                fpas.flush();
                return true;
            }
            tentativa = nova;
        }
        an->Restore();
        return false;
    };
    // avanca de t0 a t1, dividindo ao meio o passo que nao converge (ate maxCortes niveis)
    std::function<bool(REAL, REAL, int)> avanca = [&](REAL t0, REAL t1, int nivel) -> bool {
        if (passo(t0, t1, nivel)) return true;
        if (nivel >= maxCortes) return false;
        nCortes++;
        const REAL tm = 0.5 * (t0 + t1);
        return avanca(t0, tm, nivel + 1) && avanca(tm, t1, nivel + 1);
    };

    // ---- checkpoint: o estado completo (solucao, memoria de todos os pontos de integracao, face de percolacao, tempo,
    // contadores e o tamanho de cada CSV) vai para checkpoint.bin depois de cada FS e a cada CHECKPOINT_MIN minutos, com
    // troca atomica (um arquivo .tmp renomeado). Ao comecar, um checkpoint desta mesma realizacao (a assinatura abaixo)
    // e lido, os CSV voltam ao tamanho que tinham nele e o transiente continua dali; os campos nao vao no arquivo, porque
    // SEMENTE e REAL os refazem iguais. Apagado ao terminar.
    const std::string ckArq = "checkpoint.bin";
    std::string assinatura;
    {
        std::string malha = envs("MESH", "");
        malha = malha.substr(malha.find_last_of('/') + 1);
        std::ostringstream a;
        a << std::setprecision(15) << "v2 solo=" << solo << " I=" << I_mm_h << " variante=" << variante << " malha=" << malha
          << " Dt=" << Dt << " Tchuva=" << Tchuva << " Tf=" << Tf << " FS_H=";
        for (double h : fsH) a << h << ';';
        a << " fs=" << fsTol << ',' << F0 << " hl=" << hl << " hr=" << hr << " sucMax=" << sucMax << " vertical=" << verticalRain
          << " campos=" << campos << " real=" << real << " semente=" << semente << " cov=" << covC << ',' << covPhi << ','
          << covK << " L=" << Lx << ',' << Ly << " alfa=" << alfa << " M=" << Mkl << ',' << fracKL << " c=" << p.c
          << " phi=" << p.phi << " ks=" << p.ks << " neq=" << an->NEq() << " pontos=" << ips.size() << " csv=" << pwCsv
          << mecCsv << ',' << mecLimiar << ',' << bandaLimiar;
        assinatura = a.str();
    }
    namespace sfs = std::filesystem;
    auto tamanho = [](const char *a) -> int64_t {
        std::error_code ec;
        const auto s = sfs::file_size(a, ec);
        return ec ? -1 : (int64_t)s;
    };
    const int64_t magico = 0x5352464b43484b32LL; // "SRFKCHK2"
    double segAntes = 0.; // tempo das execucoes anteriores desta realizacao, ate o checkpoint de onde esta retomou
    auto salvaCheckpoint = [&](int step, REAL tt) {
        for (const auto &c : csvs) if (c.ativo) c.o->flush();
        {
            TPZBFileStream arq;
            arq.OpenWrite(ckArq + ".tmp");
            TPZStream &f = arq; // os argumentos padrao (howMany = 1) sao os de TPZStream
            f.Write(&magico);
            f.Write(&assinatura);
            const int64_t inteiros[5] = {step, nPassos, nCortes, newtonTotal, retomadas};
            f.Write(inteiros, 5);
            const double reais[3] = {tt, infAcum, segAntes + Segundos(inicio)};
            f.Write(reais, 3);
            const int64_t narq = (int64_t)csvs.size();
            f.Write(&narq);
            for (const auto &c : csvs) { const int64_t s = c.ativo ? tamanho(c.nome) : -1; f.Write(&s); }
            const int64_t ns = (int64_t)saida.size();
            f.Write(&ns);
            for (int64_t no : saida) f.Write(&no);
            const int64_t neq = an->X().Rows();
            f.Write(&neq);
            for (int64_t i = 0; i < neq; ++i) { const double v = an->X().GetVal(i, 0); f.Write(&v); }
            const int64_t nm = mem->NElements();
            f.Write(&nm);
            for (int64_t i = 0; i < nm; ++i) (*mem)[i].Write(f, 0);
            f.Write(&magico);
            arq.CloseWrite();
        }
        std::rename((ckArq + ".tmp").c_str(), ckArq.c_str());
    };
    std::vector<int64_t> tamanhosCk; // tamanho de cada CSV no checkpoint lido
    auto leCheckpoint = [&](int &step, REAL &tt) -> bool {
        if (tamanho(ckArq.c_str()) < 64) return false;
        TPZBFileStream arq;
        arq.OpenRead(ckArq);
        TPZStream &f = arq;
        int64_t m0 = 0;
        f.Read(&m0);
        std::string a;
        if (m0 == magico) f.Read(&a);
        if (m0 != magico || a != assinatura) {
            std::cerr << "checkpoint.bin de outra realizacao, parametros ou versao: ignorado\n";
            arq.CloseRead();
            return false;
        }
        int64_t inteiros[5];
        f.Read(inteiros, 5);
        double reais[3];
        f.Read(reais, 3);
        int64_t narq = 0;
        f.Read(&narq);
        std::vector<int64_t> tams(std::max<int64_t>(narq, 0));
        for (auto &s : tams) f.Read(&s);
        bool csvOk = narq == (int64_t)csvs.size();
        for (size_t k = 0; csvOk && k < csvs.size(); ++k)
            if (csvs[k].ativo && (tams[k] < 0 || tamanho(csvs[k].nome) < tams[k])) csvOk = false;
        if (!csvOk) {
            std::cerr << "checkpoint.bin sem os CSV correspondentes: ignorado\n";
            arq.CloseRead();
            return false;
        }
        int64_t ns = 0;
        f.Read(&ns);
        std::set<int64_t> s;
        for (int64_t i = 0; i < ns; ++i) { int64_t no; f.Read(&no); s.insert(no); }
        int64_t neq = 0;
        f.Read(&neq);
        if (neq != an->NEq()) { arq.CloseRead(); return false; }
        TPZFMatrix<STATE> X(neq, 1, 0.);
        for (int64_t i = 0; i < neq; ++i) { double v; f.Read(&v); X(i, 0) = v; }
        int64_t nm = 0;
        f.Read(&nm);
        if (nm != mem->NElements()) { arq.CloseRead(); return false; }
        const auto memAntes = an->SaveMemory();
        for (int64_t i = 0; i < nm; ++i) (*mem)[i].Read(f, nullptr);
        int64_t m1 = 0;
        f.Read(&m1);
        arq.CloseRead();
        if (m1 != magico) { // incompleto: volta a memoria e comeca do zero
            an->RestoreMemory(memAntes);
            std::cerr << "checkpoint.bin incompleto: ignorado\n";
            return false;
        }
        an->X() = X;
        an->Restore();
        saida = s;
        step = (int)inteiros[0];
        nPassos = inteiros[1]; nCortes = inteiros[2]; newtonTotal = inteiros[3]; retomadas = inteiros[4] + 1;
        tt = reais[0];
        infAcum = reais[1];
        segAntes = reais[2];
        tamanhosCk = tams;
        return true;
    };

    // ---- instantes: FS (e p_w em pw.csv) em FS_H; VTK_H: o ultimo estado convergido da busca do SRM (o mecanismo de
    // ruptura, com c'/F e tan(phi')/F) vai para ruptura_pontos.scal_vec.<passo>.vtk, como no driver de Rahardjo
    std::vector<REAL> Fdenso;
    for (REAL f = 0.7; f < 3.05; f += 0.1) Fdenso.push_back(std::round(f * 10.) / 10.);
    Fdenso.push_back(3.5);
    Fdenso.push_back(4.);
    const int nsteps = int(std::lround(Tf / Dt));
    std::set<int> fsSteps, vtkSteps;
    for (double h : fsH) fsSteps.insert(int(std::lround(h * hora / Dt)));
    for (double h : envhoras("VTK_H", "")) vtkSteps.insert(int(std::lround(h * hora / Dt)));

    // mecanismo do ultimo estado convergido do SRM: sqrt(J2) da deformacao plastica em cada ponto; a banda de ruptura
    // sao os pontos com sqrt(J2) >= BANDA_LIMIAR * maximo (area, centroide, extensao e propriedades medias, ponderadas
    // pela area); os pontos com sqrt(J2) >= MEC_LIMIAR * maximo vao para mecanismo.csv
    struct Banda {
        int64_t n = 0, pl = 0;
        double epsMax = 0., area = 0., xc = 0., yc = 0., xmin = 0., xmax = 0., ymin = 0., ymax = 0., c = 0., phi = 0., k = 0.;
    };
    auto analisaMecanismo = [&](REAL th, REAL F) {
        Banda b;
        std::vector<double> r(ips.size());
        for (size_t i = 0; i < ips.size(); ++i) {
            const auto &st = (*mem)[ips[i]].m_elastoplastic_state;
            TPZTensor<REAL> e(st.m_eps_p); // distorcoes de engenharia -> componentes do tensor
            e.XY() *= 0.5;
            e.XZ() *= 0.5;
            e.YZ() *= 0.5;
            r[i] = std::sqrt(std::max(REAL(0), e.J2()));
            b.epsMax = std::max(b.epsMax, r[i]);
            if (st.m_m_type != 0) b.pl++;
        }
        if (!(b.epsMax > 0.)) return b;
        b.xmin = b.ymin = 1.e300;
        b.xmax = b.ymax = -1.e300;
        for (size_t i = 0; i < ips.size(); ++i) {
            if (fmec.is_open() && r[i] >= mecLimiar * b.epsMax) {
                const auto &st = (*mem)[ips[i]].m_elastoplastic_state;
                fmec << th << ',' << F << ',' << ips[i] << ',' << xs[i] << ',' << ys[i] << ',' << r[i] << ',' << st.m_eps_p.XX()
                     << ',' << st.m_eps_p.YY() << ',' << st.m_eps_p.XY() << ',' << st.m_m_type << ',' << cPonto(i) << ','
                     << phiPonto(i) << ',' << kPonto(i) << '\n';
            }
            if (r[i] < bandaLimiar * b.epsMax) continue;
            const double a = area[i];
            b.n++;
            b.area += a;
            b.xc += a * xs[i]; b.yc += a * ys[i];
            b.c += a * cPonto(i); b.phi += a * phiPonto(i); b.k += a * kPonto(i);
            b.xmin = std::min(b.xmin, xs[i]); b.xmax = std::max(b.xmax, xs[i]);
            b.ymin = std::min(b.ymin, ys[i]); b.ymax = std::max(b.ymax, ys[i]);
        }
        b.xc /= b.area; b.yc /= b.area; b.c /= b.area; b.phi /= b.area; b.k /= b.area;
        if (fmec.is_open()) fmec.flush();
        return b;
    };

    REAL t = 0.;
    int step0 = 0;
    const bool retomou = leCheckpoint(step0, t);
    for (size_t k = 0; k < csvs.size(); ++k) {
        const auto &c = csvs[k];
        if (!c.ativo) {
            if (!retomou) std::remove(c.nome); // de uma execucao anterior com outra escolha
            continue;
        }
        if (retomou) {
            sfs::resize_file(c.nome, (uintmax_t)tamanhosCk[k]);
            c.o->open(c.nome, std::ios::app);
        } else {
            c.o->open(c.nome, std::ios::trunc);
            *c.o << c.cab << '\n';
        }
        *c.o << std::setprecision(c.prec);
    }
    if (retomou) std::cout << "retomando do checkpoint: t = " << t / hora << " h (passo " << step0 << ")" << std::endl;
    auto ultimoCk = Relogio::now();
    std::string situacao = "ok";
    for (int step = retomou ? step0 + 1 : 0; step <= nsteps; step++) {
        stepAtual = step;
        if (step > 0) {
            const REAL tn = std::min(step * Dt, Tf);
            if (!avanca(t, tn, 0)) {
                situacao = "transiente_nao_convergiu_t=" + std::to_string(t / hora);
                break;
            }
            t = tn;
        }
        if (fsSteps.count(step)) {
            const auto t0 = Relogio::now();
            const auto estado = an->X();
            const auto memoria = an->SaveMemory();
            const REAL dtAnterior = solomat->TimeStep();
            // estado hidraulico neste instante (as pressoes ficam congeladas no SRM)
            REAL pwMax = -1.e300;
            for (const auto &e : an->WaterEq()) pwMax = std::max(pwMax, estado.GetVal(e.second, 0));
            const REAL Vagua = VolumeAgua(estado);
            if (fpw.is_open()) {
                for (const auto &e : an->WaterEq()) {
                    const REAL pw = estado.GetVal(e.second, 0), pa = estado.GetVal(an->AirEq().at(e.first), 0);
                    fpw << t / hora << ',' << e.first << ',' << an->Coord(e.first, 0) << ',' << an->Coord(e.first, 1) << ','
                        << pw << ',' << pa << ',' << solomat->Hydraulics().Evaluate(pw, pa).S << ','
                        << borda.count(e.first) << ',' << saida.count(e.first) << '\n';
                }
                fpw.flush();
            }
            std::vector<TentativaSRM> tent;
            REAL Fparada = -1.;
            const REAL F = FactorSRMRetry(*an, estado, F0, fsTol, Fdenso, Fparada, tent);
            int newtonF = 0;
            for (size_t k = 0; k < tent.size(); ++k) {
                const auto &q = tent[k];
                newtonF += q.iter;
                ften << t / hora << ',' << k + 1 << ',' << q.fase << ',' << q.g << ',' << q.F << ',' << q.conv << ',' << q.iter
                     << ',' << q.retry << ',' << q.seg << '\n';
            }
            ften.flush();
            Banda b;
            REAL uMax = -1.;
            if (F > 0.) { // o ultimo estado convergido da busca, com F = FS
                b = analisaMecanismo(t / hora, F);
                uMax = 0.;
                for (const auto &e : an->DispEq())
                    uMax = std::max(uMax, std::hypot(an->NodalDisplacement(an->X(), e.first, 0),
                                                     an->NodalDisplacement(an->X(), e.first, 1)));
                if (vtkSteps.count(step)) {
                    an->SetStrengthReductionFactor(F);
                    PostProcessVTK(*an, "ruptura", step);
                }
            }
            an->SetStrengthReductionFactor(1.);
            an->SetGravityFactor(1.);
            an->SetGravityScalesSuction(true);
            an->RestoreMemory(memoria);
            an->X() = estado;
            an->Restore();
            solomat->SetTimeStep(dtAnterior);
            solomat->SetFrozenPressures(false);
            const double seg = Segundos(t0);
            auto se = [](bool ok, double v) { return ok ? Num(v) : std::string(); };
            const bool mec = F > 0., banda = b.n > 0;
            ffs << t / hora << ',' << F << ',' << Fparada << ',' << tent.size() << ',' << newtonF << ',' << se(mec, uMax) << ','
                << se(mec, b.pl) << ',' << se(mec, b.epsMax) << ',' << se(mec, b.n) << ',' << se(banda, b.area) << ','
                << se(banda, b.xc) << ',' << se(banda, b.yc) << ',' << se(banda, b.xmin) << ',' << se(banda, b.xmax) << ','
                << se(banda, b.ymin) << ',' << se(banda, b.ymax) << ',' << se(banda, b.c) << ',' << se(banda, b.phi) << ','
                << se(banda, b.k) << ',' << infAcum << ',' << Vagua << ',' << pwMax << ',' << saida.size() << ',' << seg << '\n';
            ffs.flush();
            std::cout << "t = " << t / hora << " h: FS = " << F << " (" << seg << " s)" << std::endl;
        }
        if (ckMin > 0. && step < nsteps && (fsSteps.count(step) || Segundos(ultimoCk) >= 60. * ckMin)) {
            salvaCheckpoint(step, t);
            ultimoCk = Relogio::now();
        }
    }
    for (const auto &c : csvs) if (c.ativo) c.o->close();
    std::remove(ckArq.c_str());
    const double total = segAntes + Segundos(inicio);

    // ---- resumo da realizacao: resumo.csv (cabecalho e uma linha) e resultado.txt; o FS de cada instante e a banda de
    // ruptura do FS minimo vem de fs.csv, completo tambem depois de retomadas
    std::vector<std::string> cab;
    std::vector<std::vector<std::string>> linhas;
    {
        std::ifstream f("fs.csv");
        std::string l;
        if (std::getline(f, l)) cab = SeparaCsv(l);
        while (std::getline(f, l))
            if (!l.empty()) linhas.push_back(SeparaCsv(l));
    }
    auto coluna = [&](const std::string &nome) {
        for (size_t i = 0; i < cab.size(); ++i) if (cab[i] == nome) return (int)i;
        return -1;
    };
    auto valor = [&](int lin, const std::string &nome) {
        const int c = coluna(nome);
        return (lin < 0 || c < 0 || c >= (int)linhas[lin].size()) ? std::string() : linhas[lin][c];
    };
    int iMin = -1;
    for (int i = 0; i < (int)linhas.size(); ++i)
        if (iMin < 0 || std::atof(valor(i, "FS_SRM").c_str()) < std::atof(valor(iMin, "FS_SRM").c_str())) iMin = i;
    std::vector<std::pair<std::string, std::string>> r = {
        {"real", std::to_string(real)}, {"solo", solo}, {"I_mm_h", Num(I_mm_h)}, {"variante", variante},
        {"campos", std::to_string(campos)}, {"semente", std::to_string(semente)}, {"KL_M", std::to_string(Musado)},
        {"KL_variancia", campos ? Num(fracUsada) : ""}, {"cov_c", Num(covC)}, {"cov_phi", Num(covPhi)},
        {"cov_k", Num(covK)}, {"Lx", Num(Lx)}, {"Ly", Num(Ly)}, {"alfa", Num(alfa)}, {"area_solo_m2", Num(areaSolo)},
        {"c_med", Num(ec.media)}, {"c_dp", Num(ec.dp)}, {"c_min", Num(ec.min)}, {"c_max", Num(ec.max)},
        {"phi_med", Num(ephi.media)}, {"phi_dp", Num(ephi.dp)}, {"phi_min", Num(ephi.min)}, {"phi_max", Num(ephi.max)},
        {"k_mult_med", Num(ek.media)}, {"k_mult_dp", Num(ek.dp)}, {"k_mult_min", Num(ek.min)},
        {"k_mult_max", Num(ek.max)}};
    for (double h : fsH) {
        std::string v;
        for (int i = 0; i < (int)linhas.size(); ++i)
            if (std::fabs(std::atof(valor(i, "t_h").c_str()) - h) < 1.e-6) v = valor(i, "FS_SRM");
        r.push_back({"FS_" + Num(h) + "h", v});
    }
    r.push_back({"FS_min", valor(iMin, "FS_SRM")});
    r.push_back({"t_min_h", valor(iMin, "t_h")});
    for (const char *c : {"banda_area_m2", "banda_xc", "banda_yc", "banda_ymin", "banda_c", "banda_phi_deg", "banda_k_mult"})
        r.push_back({c, valor(iMin, c)});
    int64_t newtonFS = 0;
    for (int i = 0; i < (int)linhas.size(); ++i) newtonFS += std::atol(valor(i, "newton_fs").c_str());
    r.push_back({"inf_total_m3", Num(infAcum)});
    r.push_back({"situacao", situacao});
    r.push_back({"passos", std::to_string(nPassos)});
    r.push_back({"cortes", std::to_string(nCortes)});
    r.push_back({"newton", std::to_string(newtonTotal)});
    r.push_back({"newton_fs", std::to_string(newtonFS)});
    r.push_back({"retomadas", std::to_string(retomadas)});
    r.push_back({"seg", Num(total)});
    r.push_back({"memoria_MB", Num(PicoMemoriaMB())});
    {
        std::ofstream o("resumo.csv");
        for (size_t i = 0; i < r.size(); ++i) o << (i ? "," : "") << r[i].first;
        o << '\n';
        for (size_t i = 0; i < r.size(); ++i) o << (i ? "," : "") << r[i].second;
        o << '\n';
    }
    const std::string fsMin = valor(iMin, "FS_SRM"), tMin = valor(iMin, "t_h");
    std::ofstream res("resultado.txt");
    res << "real=" << real << " FS_min=" << fsMin << " t_min_h=" << tMin << " situacao=" << situacao << " passos=" << nPassos
        << " cortes=" << nCortes << " newton=" << newtonTotal << " seg=" << Num(total) << '\n';
    std::cout << "real " << real << ": FS_min = " << fsMin << " em t = " << tMin << " h, " << situacao << ", " << total << " s"
              << std::endl;
    delete an;
    DeleteCompMesh(cmesh);
    delete gmesh;
    return situacao == "ok" ? 0 : 1;
}
