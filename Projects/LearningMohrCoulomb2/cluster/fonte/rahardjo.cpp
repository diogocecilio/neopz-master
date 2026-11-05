// rahardjo.cpp - Rahardjo et al. (2007), Serie A, alpha = 45 graus, rodado com o codigo de LearningMohrCoulomb2 SEM alterar o projeto.
// Inclui o main.cpp do projeto (FactorSRM, FactorGIM, PostProcessVTK, CreateSlopeCompMesh...) e replica o laco da
// chuva de main() linha a linha, acrescentando: parametros por variavel de ambiente, balanco de massa por passo,
// dump do estado (X e memoria) e casos de FS isolados.
#define main projeto_main_original
#include "main.cpp"
#undef main

#include <chrono>
#include <cstdlib>
#include <sstream>

static double envd(const char *k, double d) { const char *v = std::getenv(k); return v ? std::atof(v) : d; }
static int envi(const char *k, int d) { const char *v = std::getenv(k); return v ? std::atoi(v) : d; }
static std::string envs(const char *k, const std::string &d) { const char *v = std::getenv(k); return v ? std::string(v) : d; }
static std::set<int> envlist(const char *k) {
    std::set<int> s; const char *v = std::getenv(k); if (!v) return s;
    std::stringstream ss(v); std::string t; while (std::getline(ss, t, ',')) if (!t.empty()) s.insert(std::atoi(t.c_str()));
    return s;
}

static void WriteTensor(std::ofstream &o, const TPZTensor<REAL> &t) {
    for (int i = 0; i < 6; i++) { const double v = t.fData[i]; o.write((const char *)&v, sizeof(double)); }
}
// Dump binario: neq, X[neq]
static void DumpState(TPZUnsatAnalysis &an, TUnsat *, int step) {
    std::ofstream o("state_" + std::to_string(step) + ".bin", std::ios::binary);
    const int64_t neq = an.NEq();
    o.write((const char *)&neq, sizeof(int64_t));
    for (int64_t i = 0; i < neq; i++) { const double v = an.X().GetVal(i, 0); o.write((const char *)&v, sizeof(double)); }
}

static double RssMB() {
    std::ifstream f("/proc/self/statm"); long a = 0, b = 0; f >> a >> b; return double(b) * 4096. / 1048576.;
}

// FactorGIM do projeto + nova tentativa do lambda que falhou, com incremento < tol
static REAL FactorGIMRetry(TPZUnsatAnalysis &an, const TPZFMatrix<STATE> &x, REAL tol, const std::vector<REAL> &lams, REAL &Lparada)
{
    an.ResetMemory();
    TPZFMatrix<STATE> &X = an.X();
    X = x;
    for (int64_t i = 0; i < an.NEq(); i++)
        if (!an.IsPressure(i)) X(i, 0) = 0.;
    an.Restore();
    const auto fixed = an.FrozenFixed(X);
    an.SetStrengthReductionFactor(1.);
    an.SetGravityScalesSuction(false);
    REAL lok = 0., lfail = -1.;
    Lparada = -1.;
    size_t next = 0;
    while (true) {
        REAL lam;
        bool retry = false;
        if (lfail > 0.) {
            if (lfail - lok < tol) {
                if (Lparada < 0.) Lparada = lok;
                lam = lfail; retry = true;
            } else lam = 0.5 * (lok + lfail);
        } else {
            while (next < lams.size() && lams[next] <= lok) next++;
            if (next == lams.size()) break;
            lam = lams[next];
        }
        an.SetGravityFactor(lam);
        const auto r = an.Solve(0., fixed, {}, true);
        if (r.conv) {
            an.Accept();
            lok = lam;
            if (retry) lfail = -1.;
        } else {
            an.Restore();
            if (retry) break;
            lfail = lam;
        }
    }
    if (Lparada < 0.) Lparada = lok;
    an.SetGravityScalesSuction(true);
    an.SetGravityFactor(1.);
    return lok;
}

// FactorSRM do projeto + nova tentativa do F que falhou, agora com incremento < tol; Fparada = valor em que o projeto para
static REAL FactorSRMRetry(TPZUnsatAnalysis &an, const TPZFMatrix<STATE> &x, REAL F0, REAL tol, const std::vector<REAL> &Fs, REAL &Fparada)
{
    an.ResetMemory();
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
        if (!an.Solve(0., fixed, {}, true).conv) { Fparada = -1.; return -1.; }
        an.Accept();
    }
    REAL Fok = F0, Ffail = -1.;
    Fparada = -1.;
    size_t next = 0;
    while (true) {
        REAL F;
        bool retry = false;
        if (Ffail > 0.) {
            if (Ffail - Fok < tol) {
                if (Fparada < 0.) Fparada = Fok;
                F = Ffail; retry = true;
            } else F = 0.5 * (Fok + Ffail);
        } else {
            while (next < Fs.size() && Fs[next] <= Fok) next++;
            if (next == Fs.size()) break;
            F = Fs[next];
        }
        an.SetStrengthReductionFactor(F);
        const auto r = an.Solve(0., fixed, {}, true);
        std::cout << "SRMr: F = " << F << (r.conv ? " convergiu" : " nao convergiu") << " (" << r.iter << " iteracoes)" << (retry ? " [nova tentativa]" : "") << std::endl;
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
    const int porder = 2;
    const REAL hora = 3600., dia = 24. * hora;
    const REAL I_mm_h = envd("I_MM_H", 36.), Ir = I_mm_h / (1000. * hora);
    const REAL Tchuva = envd("TCHUVA_H", 24.) * hora, Tf = envd("TF_H", 240.) * hora, Dt = envd("DT_H", 0.1) * hora;
    const REAL xl = envd("XL", 0.), xr = envd("XR", 70.);
    const REAL Hw = envd("HW", 5.), beta = envd("BETA_DEG", 7.) * M_PI / 180., sucMax = envd("SUCMAX", 75.);
    const REAL xtoe = envd("XTOE", 40.), ytoe = envd("YTOE", 20.);
    const REAL hl = envd("HL", ytoe - Hw + (xtoe - xl) * std::tan(beta)), hr = envd("HR", ytoe - Hw + (xtoe - xr) * std::tan(beta));
    const int cadaFS = envi("CADA_FS", 10); // 0: nunca
    std::set<int> fsSteps = envlist("FS_STEPS"), dumpSteps = envlist("DUMP_STEPS"), vtkSteps, gimSteps;
    const int comGIM = envi("GIM", 0);
    const int vtk = envi("VTK", 0), verticalRain = envi("CHUVA_VERTICAL", 0);
    const REAL fsTol = envd("FS_TOL", 0.01), F0 = envd("FS_F0", 0.6);
    const REAL inicioDetalhe = 0., fimDetalhe = Tchuva, tolTempo = 1.e-8 * hora;

    TPZGeoMesh *gmeshslope = ReadGmshMesh(envs("MESH", "malhas/rahardjo_a45.msh"), {{"Solo", matSolo}},
        {{"Base", bcBase}, {"Direita", bcDireita}, {"Superficie", bcSuperficie}, {"Esquerda", bcEsquerda}});
    if (vtk) { std::ofstream out2("talude.vtk"); TPZVTKGeoMesh::PrintGMeshVTK(gmeshslope, out2, true); }
    // Tabela 1 e resistencia do artigo; E, nu, ar e agua como na Tabela de parametros do documento (hipoteses)
    const std::string solo = envs("SOLO", "f50");
    const REAL aFX = solo == "f10" ? 10. : solo == "f50" ? 50. : solo == "f100" ? 100. : -1.;
    if (aFX < 0.) { std::cerr << "SOLO deve ser f10, f50 ou f100\n"; return 2; }
    const std::string variante = envs("VARIANTE", "iii"); // i: completo; ii: ar passivo; iii: desacoplado, ar passivo, esqueleto rigido
    UnsatParams p;
    p.E = 2.e4, p.nu = 0.3, p.c = 10., p.phi = 26. * M_PI / 180., p.gamma = 20.;
    p.n = 0.45, p.ks = solo == "f10" ? 1.e-4 : solo == "f50" ? 1.e-5 : 1.e-6;
    p.hydro.fRetention = TPZUnsaturatedHydraulics::EFredlundXing;
    p.hydro.fKrw = TPZUnsaturatedHydraulics::EPowerLaw;
    p.hydro.fKra = TPZUnsaturatedHydraulics::EBrooksCorey;
    p.hydro.fA = aFX, p.hydro.fN = 1., p.hydro.fM = 1., p.hydro.fPower = 4.;
    p.bishop = TUnsat::EChiConstant, p.chi = 1.; // phi^b = phi'
    p.coupled = variante == "i" || variante == "ii";
    p.E = envd("E", p.E); p.nu = envd("NU", p.nu); p.c = envd("COES", p.c); p.phi = envd("PHI_DEG", p.phi * 180. / M_PI) * M_PI / 180.;
    p.gamma = envd("GAMMA", p.gamma); p.ks = envd("KS", p.ks); p.n = envd("POROS", p.n); p.Kw = envd("KW", p.Kw);
    p.mixture = envi("MIXTURE", p.mixture); p.lumped = envi("LUMPED", p.lumped);
    TPZMultiphysicsCompMesh *cslope = CreateSlopeCompMesh(gmeshslope, porder, p);
    // solver linear: PARDISO (padrao, com MKL) ou SOLVER=skyline, o LU usado por todas as rodadas anteriores
    TPZUnsatAnalysis::SetUseSkyline(envs("SOLVER", "pardiso") == "skyline");
    auto *anslope = new TPZUnsatAnalysis(cslope, {{bcBase, {true, true}}, {bcDireita, {true, false}}, {bcEsquerda, {true, false}}});
    // pardiso-metis (padrao) ou pardiso-nativa (PZ_PARDISO_ORDEM=nativa, a ordem das rodadas anteriores, registradas so como "pardiso")
    const char *ordemPardiso = std::getenv("PZ_PARDISO_ORDEM");
    const std::string solverTxt = !TPZUnsatAnalysis::UsesPardiso() ? "skyline"
        : (ordemPardiso && std::string(ordemPardiso) == "nativa") ? "pardiso-nativa" : "pardiso-metis";
    std::cout << "solver: " << solverTxt << std::endl;
    anslope->SetTolerances(envd("TOLU", 1.e-7), envd("TOLP", 1.e-6), envi("MAXIT", 40), envi("MAXIT_FROZEN", 25));
    auto *reservatorio = dynamic_cast<TUnsat *>(cslope->FindMaterial(matSolo));
    if (!reservatorio) DebugStop();
    const auto borda = anslope->BoundaryNodes(bcSuperficie, verticalRain != 0);
    if (borda.empty()) DebugStop();

    // Geometria auxiliar: nos, triangulos e massa concentrada de cada vertice (area/3)
    std::map<int64_t, REAL> massa;
    {
        std::ofstream fn("nodes.csv"), ft("tri.csv");
        fn << "node,x,y,eqW,eqA,eqU\n" << std::setprecision(15);
        for (const auto &e : anslope->WaterEq())
            fn << e.first << ',' << anslope->Coord(e.first, 0) << ',' << anslope->Coord(e.first, 1) << ',' << e.second << ','
               << anslope->AirEq().at(e.first) << ',' << anslope->DispEq().at(e.first) << '\n';
        ft << "n0,n1,n2\n";
        for (int64_t i = 0; i < gmeshslope->NElements(); i++) {
            TPZGeoEl *g = gmeshslope->Element(i);
            if (!g || g->Dimension() != 2) continue;
            const REAL x0 = anslope->Coord(g->NodeIndex(0), 0), y0 = anslope->Coord(g->NodeIndex(0), 1);
            const REAL x1 = anslope->Coord(g->NodeIndex(1), 0), y1 = anslope->Coord(g->NodeIndex(1), 1);
            const REAL x2 = anslope->Coord(g->NodeIndex(2), 0), y2 = anslope->Coord(g->NodeIndex(2), 1);
            const REAL A = 0.5 * std::fabs((x1 - x0) * (y2 - y0) - (x2 - x0) * (y1 - y0));
            for (int k = 0; k < 3; k++) massa[g->NodeIndex(k)] += A / 3.;
            ft << g->NodeIndex(0) << ',' << g->NodeIndex(1) << ',' << g->NodeIndex(2) << '\n';
        }
    }
    auto VolumeAgua = [&](const TPZFMatrix<STATE> &x) {
        REAL v = 0.;
        for (const auto &m : massa) {
            const REAL pw = x.GetVal(anslope->WaterEq().at(m.first), 0), pa = x.GetVal(anslope->AirEq().at(m.first), 0);
            v += m.second * p.n * reservatorio->Hydraulics().Evaluate(pw, pa).S;
        }
        return v;
    };

    // ---- modo FS: le estados gravados pelo transiente e calcula os FS (so dependem das pressoes), em ordem de prioridade
    if (const char *lista = std::getenv("FS_ESTADOS")) {
        const std::string nomeCsv = envs("FS_CSV", "fs.csv");
        auto chave = [](double h) { return std::lround(h * 1000.); };
        // retomavel: instantes ja gravados em fs.csv sao pulados
        std::set<long> feitos, gimH, vtkH, projetoH;
        bool novo = true;
        { std::ifstream ant(nomeCsv); std::string l; if (std::getline(ant, l)) { novo = false; while (std::getline(ant, l)) if (!l.empty()) feitos.insert(chave(std::atof(l.c_str()))); } }
        std::ofstream fc(nomeCsv, std::ios::app);
        if (novo) fc << "t_h,FS_SRM,FS_SRM_1a_parada,FS_projeto,FS_GIM,FS_GIM_1a_parada,seg\n";
        fc << std::setprecision(10);
        auto leH = [&](const char *k, std::set<long> &alvo) { const char *v = std::getenv(k); if (!v) return; std::stringstream ss(v); std::string tk; while (std::getline(ss, tk, ',')) if (!tk.empty()) alvo.insert(chave(std::atof(tk.c_str()))); };
        leH("GIM_H", gimH); leH("VTK_H", vtkH); leH("FS_PROJETO_H", projetoH);
        std::vector<REAL> Fdenso; for (REAL f = 0.7; f < 3.05; f += 0.1) Fdenso.push_back(std::round(f * 10.) / 10.); Fdenso.push_back(3.5); Fdenso.push_back(4.);
        const std::vector<REAL> lams = {0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 6.0, 7.0, 8.0, 10.0};
        std::stringstream ss(lista); std::string tk;
        while (std::getline(ss, tk, ',')) {
            if (tk.empty()) continue;
            const double h = std::atof(tk.c_str());
            if (feitos.count(chave(h))) continue;
            const int passo = int(std::lround(h * hora / Dt));
            std::ifstream e("state_" + std::to_string(passo) + ".bin", std::ios::binary);
            if (!e) { std::cerr << "sem estado para t = " << h << " h\n"; continue; }
            int64_t neq = 0; e.read((char *)&neq, sizeof(int64_t));
            if (neq != cslope->NEquations()) { std::cerr << "neq diferente\n"; return 3; }
            TPZFMatrix<STATE> xe(neq, 1, 0.);
            for (int64_t i = 0; i < neq; i++) { double v; e.read((char *)&v, sizeof(double)); xe(i, 0) = v; }
            const auto t0 = std::chrono::steady_clock::now();
            REAL Fp, Lp = -1., G = -1., Fproj = -1.;
            const REAL F = FactorSRMRetry(*anslope, xe, F0, fsTol, Fdenso, Fp);
            if (vtkH.count(chave(h)) && F > 0.) { anslope->SetStrengthReductionFactor(F); PostProcessVTK(*anslope, "mec", passo); anslope->SetStrengthReductionFactor(1.); }
            if (projetoH.count(chave(h))) Fproj = FactorSRM(*anslope, xe, F0, fsTol);
            if (gimH.count(chave(h))) G = FactorGIMRetry(*anslope, xe, fsTol, lams, Lp);
            const double seg = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
            fc << h << ',' << F << ',' << Fp << ',' << Fproj << ',' << G << ',' << Lp << ',' << seg << '\n';
            fc.flush();
        }
        return 0;
    }

    // ---- estado inicial (identico a main(); SUC_UNIF/SECO sao os casos isolados de FS)
    TPZFMatrix<STATE> xslope(cslope->NEquations(), 1, 0.);
    const REAL sucUnif = envd("SUC_UNIF", -1.);
    const int seco = envi("SECO", 0);
    for (const auto &eq : anslope->WaterEq()) {
        const REAL x = anslope->Coord(eq.first, 0), y = anslope->Coord(eq.first, 1);
        const REAL hf = hl + (hr - hl) * (x - xl) / (xr - xl);
        xslope(eq.second, 0) = std::max(p.gw * (hf - y), -sucMax);
        if (sucUnif >= 0.) xslope(eq.second, 0) = -sucUnif;
        if (seco) xslope(eq.second, 0) = 0.;
    }
    if (envi("FS_ONLY", 0)) { // FS de um estado imposto, sem transiente
        if (const char *arq = std::getenv("PW_CSV")) {
            std::ifstream f(arq);
            REAL dmax = 0.; int lidos = 0;
            std::string lin;
            while (std::getline(f, lin)) {
                double a, b, c;
                if (std::sscanf(lin.c_str(), "%lf,%lf,%lf", &a, &b, &c) != 3) continue;
                int64_t melhor = -1; REAL dm = 1.e30;
                for (const auto &e : anslope->WaterEq()) {
                    const REAL d = std::hypot(anslope->Coord(e.first, 0) - a, anslope->Coord(e.first, 1) - b);
                    if (d < dm) { dm = d; melhor = e.first; }
                }
                xslope(anslope->WaterEq().at(melhor), 0) = c;
                dmax = std::max(dmax, dm); lidos++;
            }
            std::cout << "PW_CSV: " << lidos << " nos, maior distancia " << dmax << std::endl;
        }
        anslope->X() = xslope;
        anslope->Restore();
        if (envi("FS_RETRY", 0)) {
            REAL Fp;
            const REAL Fr = FactorSRMRetry(*anslope, xslope, F0, fsTol, {0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.5, 3.0, 3.5, 4.0}, Fp);
            std::cout << "RESULTADO RETRY FS = " << std::setprecision(10) << Fr << " parada do projeto = " << Fp << std::endl;
            return 0;
        }
        const REAL F = FactorSRM(*anslope, xslope, F0, fsTol);
        const REAL L = FactorGIM(*anslope, xslope, fsTol);
        std::cout << "RESULTADO FS_SRM " << std::setprecision(10) << F << " FS_GIM " << L << std::endl;
        return 0;
    }
    anslope->X() = xslope;
    anslope->Restore();
    const auto inicial = anslope->FrozenFixed(xslope);
    anslope->SetStrengthReductionFactor(1.);
    anslope->SetGravityScalesSuction(true);
    bool convergiu = true;
    for (int k = 1; k <= 4; k++) {
        anslope->SetGravityFactor(REAL(k) / 4.);
        if (!anslope->Solve(0., inicial, {}, true).conv) { convergiu = false; break; }
        anslope->Accept();
    }
    if (convergiu) anslope->ResetDisplacements();
    else std::cerr << "Estado inicial nao convergiu com a resistencia real.\n";

    // variante iii: esqueleto rigido, so a percolacao e resolvida (SEEP/W); as variantes i e ii resolvem o equilibrio
    auto contorno = variante == "iii" ? anslope->AllDisplacements(anslope->X()) : anslope->DisplacementConstraints();
    std::set<int64_t> eqEsq, eqDir;
    for (const int id : {bcEsquerda, bcDireita}) {
        const REAL h = id == bcEsquerda ? hl : hr;
        for (const auto &b : anslope->BoundaryNodes(id)) {
            const REAL y = anslope->Coord(b.first, 1);
            if (y <= h) {
                contorno[anslope->WaterEq().at(b.first)] = p.gw * (h - y);
                (id == bcEsquerda ? eqEsq : eqDir).insert(anslope->WaterEq().at(b.first));
            }
        }
    }
    for (const auto &b : borda) contorno[anslope->AirEq().at(b.first)] = 0.;
    if (variante != "i") for (const auto &e : anslope->AirEq()) contorno[e.second] = 0.; // ar passivo: p_a = 0 em todo o dominio
    {
        std::ofstream c("caso.txt");
        c << "solo=" << solo << " variante=" << variante << " I_mm_h=" << I_mm_h << " Tchuva_h=" << Tchuva / hora << " Tf_h=" << Tf / hora
          << " Dt_h=" << Dt / hora << " hl=" << hl << " hr=" << hr << " sucMax=" << sucMax << " neq=" << anslope->NEq()
          << " chuva_vertical=" << verticalRain << " E=" << p.E << " nu=" << p.nu << " c=" << p.c << " phi=" << p.phi * 180. / M_PI
          << " a=" << p.hydro.fA << " ks=" << p.ks << " solver=" << solverTxt << '\n';
    }

    std::ofstream csv("chuva.csv");
    csv << "passo,tempo_s,tempo_h,tempo_dias,chuva_media_passo_mm_h,FS_SRM,FS_GIM,seg_FS,rss_MB\n" << std::setprecision(12);
    std::ofstream bal("balanco.csv");
    bal << "passo,t_h,chuva_mm_h,chuva_pot,inf_fluxo,inf_sat,esq,dir,dV_lump,dV_acop,dV_comp,somaR,residuo_livre,nsat,"
           "its_ativo,newton,normp,eu,pl_pts,pw_min,pw_max,pa_min,pa_max,seg_passo\n" << std::setprecision(15);
    std::set<int64_t> saida;
    const REAL tolP = 1.e-6, tolQ = 1.e-10;
    // tempos: Dt ate T2, depois Dt2 (T2 = Tf reproduz o passo uniforme de main())
    const REAL T2 = std::min(envd("T2_H", Tf / hora) * hora, Tf), Dt2 = envd("DT2_H", Dt / hora) * hora;
    std::vector<REAL> tempos;
    {
        const int n1 = int(std::lround(T2 / Dt));
        for (int i = 0; i <= n1; i++) tempos.push_back(i * Dt);
        const int n2 = int(std::ceil((Tf - tempos.back()) / Dt2 - 1.e-9));
        const REAL T1 = tempos.back();
        for (int j = 1; j <= n2; j++) tempos.push_back(std::min(T1 + j * Dt2, Tf));
    }
    auto horas = [&](const char *k, std::set<int> &alvo) { // instantes pedidos em horas -> passo mais proximo
        const char *v = std::getenv(k);
        if (!v) return;
        std::stringstream ss(v);
        std::string tk;
        while (std::getline(ss, tk, ',')) {
            if (tk.empty()) continue;
            const REAL h = std::atof(tk.c_str()) * hora;
            int melhor = 0;
            for (int i = 0; i < int(tempos.size()); i++)
                if (std::fabs(tempos[i] - h) < std::fabs(tempos[melhor] - h)) melhor = i;
            if (std::fabs(tempos[melhor] - h) < 1.e-3 * hora) alvo.insert(melhor);
        }
    };
    horas("FS_H", fsSteps); horas("DUMP_H", dumpSteps); horas("VTK_H", vtkSteps); horas("GIM_H", gimSteps);
    const int nsteps = int(tempos.size()) - 1;
    REAL t = 0., chuva = 0.;
    if (dumpSteps.count(0)) DumpState(*anslope, reservatorio, 0);
    for (int step = 0; step <= nsteps && convergiu; step++) {
        if (step > 0) {
            const auto t0 = std::chrono::steady_clock::now();
            const REAL tn = tempos[step], dt = tn - t;
            const REAL chuvaPasso = Ir * std::max(REAL(0), std::min(tn, Tchuva) - t) / dt;
            auto tentativa = saida;
            for (const auto &b : borda) {
                const int64_t node = b.first;
                if (anslope->X().GetVal(anslope->WaterEq().at(node), 0) > tolP) tentativa.insert(node);
            }
            convergiu = false;
            TPZUnsatAnalysis::StepResult sr;
            int itsAtivo = 0;
            TPZFMatrix<STATE> Rfinal;
            std::map<int64_t, REAL> fqFinal;
            for (int it = 0; it < 40; it++) {
                itsAtivo = it + 1;
                auto fixed = contorno;
                std::map<int64_t, REAL> fq;
                for (const auto &b : borda) {
                    const int64_t node = b.first, eq = anslope->WaterEq().at(node);
                    if (tentativa.count(node)) fixed[eq] = 0.;
                    else fq[eq] = chuvaPasso * b.second;
                }
                sr = anslope->Solve(dt, fixed, fq, false);
                if (!sr.conv) break;
                const auto R = anslope->FullResidual();
                const TPZFMatrix<STATE> solucao = cslope->Solution();
                auto nova = tentativa;
                for (const auto &b : borda) {
                    const int64_t node = b.first, eq = anslope->WaterEq().at(node);
                    if (tentativa.count(node)) {
                        if (R.GetVal(eq, 0) / (dt * b.second) > chuvaPasso + tolQ) nova.erase(node);
                    } else if (solucao.GetVal(eq, 0) > tolP) nova.insert(node);
                }
                if (nova == tentativa) { convergiu = true; Rfinal = R; fqFinal = fq; break; }
                tentativa = nova;
            }
            if (!convergiu) {
                anslope->Restore();
                std::cerr << "Passo " << step << " nao convergiu (t = " << t / hora << " h; its_ativo=" << itsAtivo
                          << ", newton=" << sr.iter << ", eu=" << sr.eu << ", normp=" << sr.normp << ")\n";
                break;
            }
            // ---- balanco de massa do passo (antes de Accept, com a solucao convergida carregada na malha)
            const TPZFMatrix<STATE> xnovo = cslope->Solution();
            const REAL V0 = VolumeAgua(anslope->X()), V1 = VolumeAgua(xnovo);
            reservatorio->SetHydroMechanicalCoupling(false);
            const auto Rsem = anslope->FullResidual();
            reservatorio->SetHydroMechanicalCoupling(p.coupled);
            reservatorio->SetWaterBulkModulus(1.e40);
            const auto Rinc = anslope->FullResidual();
            reservatorio->SetWaterBulkModulus(p.Kw);
            REAL somaR = 0., somaSem = 0., somaInc = 0.;
            for (const auto &e : anslope->WaterEq()) {
                somaR += Rfinal.GetVal(e.second, 0);
                somaSem += Rsem.GetVal(e.second, 0);
                somaInc += Rinc.GetVal(e.second, 0);
            }
            REAL chuvaPot = 0., infFluxo = 0., infSat = 0.;
            for (const auto &b : borda) {
                const int64_t eq = anslope->WaterEq().at(b.first);
                chuvaPot += chuvaPasso * b.second * dt;
                if (tentativa.count(b.first)) infSat += Rfinal.GetVal(eq, 0);
                else infFluxo += chuvaPasso * b.second * dt;
            }
            REAL qEsq = 0., qDir = 0.;
            for (auto e : eqEsq) qEsq += Rfinal.GetVal(e, 0);
            for (auto e : eqDir) qDir += Rfinal.GetVal(e, 0);
            REAL pwmin = 1e300, pwmax = -1e300, pamin = 1e300, pamax = -1e300, residuoLivre = 0.;
            for (const auto &e : anslope->WaterEq()) {
                const REAL pw = xnovo.GetVal(e.second, 0), pa = xnovo.GetVal(anslope->AirEq().at(e.first), 0);
                pwmin = std::min(pwmin, pw); pwmax = std::max(pwmax, pw);
                pamin = std::min(pamin, pa); pamax = std::max(pamax, pa);
            }
            residuoLivre = somaR - infFluxo - infSat - qEsq - qDir;
            anslope->Accept();
            saida = tentativa;
            t = tn;
            chuva = chuvaPasso;
            const double seg = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
            bal << step << ',' << t / hora << ',' << chuva * 1000. * hora << ',' << chuvaPot << ',' << infFluxo << ',' << infSat
                << ',' << qEsq << ',' << qDir << ',' << V1 - V0 << ',' << somaR - somaSem << ',' << somaR - somaInc << ','
                << somaR << ',' << residuoLivre << ',' << saida.size() << ',' << itsAtivo << ',' << sr.iter << ',' << sr.normp
                << ',' << sr.eu << ',' << reservatorio->NPlasticPoints() << ',' << pwmin << ',' << pwmax << ',' << pamin << ','
                << pamax << ',' << seg << '\n';
            bal.flush();
            std::cout << "Passo " << step << ": t = " << t / hora << " h; chuva media = " << chuva * 1000. * hora << " mm/h\n";
        }
        if (dumpSteps.count(step) && step > 0) DumpState(*anslope, reservatorio, step);

        const bool detalhar = t >= inicioDetalhe - tolTempo && t <= fimDetalhe + tolTempo;
        if ((cadaFS > 0 && step % cadaFS == 0) || fsSteps.count(step) || (cadaFS > 0 && step == nsteps)) {
            const auto tf0 = std::chrono::steady_clock::now();
            const auto estado = anslope->X();
            const auto memoria = anslope->SaveMemory();
            const REAL dtAnterior = reservatorio->TimeStep();
            const REAL FSslope = FactorSRM(*anslope, estado, F0, fsTol);
            if (vtkSteps.count(step) && FSslope > 0.) {
                anslope->SetStrengthReductionFactor(FSslope);
                PostProcessVTK(*anslope, "talude_chuva_srm", step);
                anslope->SetStrengthReductionFactor(1.);
            }
            const REAL FGslope = (comGIM || gimSteps.count(step)) ? FactorGIM(*anslope, estado, fsTol) : REAL(-1.);
            if (vtkSteps.count(step) && FGslope > 0.) {
                anslope->SetGravityScalesSuction(false);
                anslope->SetGravityFactor(FGslope);
                PostProcessVTK(*anslope, "talude_chuva_gim", step);
            }
            anslope->SetStrengthReductionFactor(1.);
            anslope->SetGravityFactor(1.);
            anslope->SetGravityScalesSuction(true);
            anslope->RestoreMemory(memoria);
            anslope->X() = estado;
            anslope->Restore();
            reservatorio->SetTimeStep(dtAnterior);
            reservatorio->SetFrozenPressures(false);
            const double segFS = std::chrono::duration<double>(std::chrono::steady_clock::now() - tf0).count();
            csv << step << ',' << t << ',' << t / hora << ',' << t / dia << ',' << chuva * 1000. * hora << ',' << FSslope << ','
                << FGslope << ',' << segFS << ',' << RssMB() << '\n';
            csv.flush();
            std::cout << "  FS (SRM) = " << FSslope << "; FS (GIM) = " << FGslope << '\n';
            if (vtkSteps.count(step)) PostProcessVTK(*anslope, "talude_chuva", step);
        }
    }
    delete anslope;
    DeleteCompMesh(cslope);
    delete gmeshslope;
    return convergiu ? 0 : 1;
}
