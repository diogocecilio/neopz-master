// V5_EncostaDuasCamadas: testes Catch2 (marca [V5]; os demorados tambem [.lento], fora do ctest padrao). Referencia:
// Abbasov et al. (2024), Comput. Geotech. 167:106102 (dados/ e README.md). Os testes padrao usam a malha grossa
// (malhas/encosta_grossa.msh) e um estado de 10 anos calculado depressa (dt de 1 dia dobrando ate 1 ano: o permanente,
// que o artigo chama de "equilibrium"/"steady state", e alcancado em ~3 anos); o [.lento] roda a malha do caso completo
// com os passos do main.cpp.
#define CATCH_CONFIG_MAIN
#include <catch2/catch.hpp>

#include "problema.h"

#include <map>
#include <memory>

using namespace v5;

namespace {

/// valores digitalizados de um mapa (dados/figN_*.csv): t_h, x_zona, z_zona, valor
struct PontoMapa { double t, x, z, v; };
std::vector<PontoMapa> LeMapa(const std::string &arq) {
    std::ifstream f(std::string(VERIF_DIR) + "/dados/" + arq);
    std::vector<PontoMapa> r;
    std::string l;
    while (std::getline(f, l)) {
        if (l.empty() || l[0] == '#' || l[0] == 't') continue;
        PontoMapa p;
        if (std::sscanf(l.c_str(), "%lf,%lf,%lf,%lf", &p.t, &p.x, &p.z, &p.v) == 4) r.push_back(p);
    }
    return r;
}

/// comparacao de um mapa num instante: diferenca media e RMS nos pontos digitalizados (o nosso valor amostrado no
/// mesmo ponto); so os pontos a mais de 'margem' da interface (a mancha de cor da linha preta e a transicao brusca)
struct DifMapa { double media = 0., rms = 0.; int n = 0; };

/// uma execucao (etapa 1 e etapa 2) com o que os testes usam de cada instante
struct Rodada {
    Config cf;
    std::unique_ptr<Encosta> e;
    bool ok1 = false, ok2 = false;
    double seg = 0.;
    InfoPasso fim1;
    std::map<double, Instante> inst;
    std::map<double, DifMapa> difTheta, difSigs, difLFS;
};

DifMapa Compara(Encosta &e, const std::vector<PontoMapa> &mapa, double t, int var, double margem = 0.15) {
    DifMapa d;
    double s = 0., s2 = 0.;
    for (const auto &p : mapa) {
        if (std::fabs(p.t - t) > 1.e-9) continue;
        const double x = p.x + kXZona0, z = p.z + kZZona0;
        if (std::fabs(z - CotaInterface(x)) < margem) continue;
        const auto a = e.Amostrar(x, z);
        if (a.el < 0) continue;
        const double v = var == 0 ? a.theta : var == 1 ? a.sigs : std::min(a.lfs, 3.);
        s += v - p.v;
        s2 += (v - p.v) * (v - p.v);
        d.n++;
    }
    if (d.n) d.media = s / d.n, d.rms = std::sqrt(s2 / d.n);
    return d;
}

/// roda um caso: rapido = malha grossa, etapa 1 com dt de 1 dia dobrando ate 1 ano e etapa 2 com dt = 0,25 h;
/// senao a malha e os passos do caso completo (main.cpp)
Rodada &Roda(const std::string &caso, bool rapido) {
    static std::map<std::string, std::unique_ptr<Rodada>> cache;
    const std::string chave = caso + (rapido ? "_rapido" : "_completo");
    auto it = cache.find(chave);
    if (it != cache.end()) return *it->second;
    auto r = std::make_unique<Rodada>();
    const auto t0 = Relogio::now();
    r->cf = Config::CasoPadrao(caso);
    if (rapido) {
        r->cf.malha = "encosta_grossa.msh";
        r->cf.dt0_h = 24., r->cf.fatorDt = 2., r->cf.dtMax_d = 365.;
        r->cf.dt2_h = 0.25;
    }
    Verboso() = 0;
    r->e = std::make_unique<Encosta>(r->cf, LeMalha(CaminhoMalha(r->cf.malha)));
    Encosta &e = *r->e;
    e.CondicaoInicial();
    r->ok1 = e.Etapa1();
    r->fim1 = e.ultimo;
    static const auto th4 = LeMapa("fig4_theta.csv"), sg5 = LeMapa("fig5_sigs.csv"), lf6 = LeMapa("fig6_lfs.csv");
    if (r->ok1)
        r->ok2 = e.Etapa2([&](double th) {
            if (!e.Tensoes()) return;
            r->inst[th] = AvaliaInstante(e, th);
            r->difTheta[th] = Compara(e, th4, th, 0);
            r->difSigs[th] = Compara(e, sg5, th, 1);
            r->difLFS[th] = Compara(e, lf6, th, 2);
        });
    r->seg = Segundos(t0);
    std::cout << "V5 " << chave << ": " << e.passos << " passos, " << r->seg << " s" << std::endl;
    for (const auto &kv : r->inst) {
        const Instante &i = kv.second;
        std::cout << "  t = " << kv.first << " h: LFS medio " << i.m.lfsMedio << " (limitado a 3: " << i.m.lfsMedioLim3
                  << "), A(LFS < 1) " << i.m.aLFS << " m2, nivel na crista " << i.nivelCrista << " m, frente "
                  << i.frenteCrista << " m, theta sup crista " << i.thetaSupCrista << ", succao sup crista "
                  << i.sigsSupCrista << " kPa; mapas (RMS): theta " << r->difTheta[kv.first].rms << ", succao "
                  << r->difSigs[kv.first].rms << " kPa, LFS (lim. 3) " << r->difLFS[kv.first].rms << std::endl;
    }
    return *cache.emplace(chave, std::move(r)).first->second;
}

// valores de referencia do artigo (Fig. 14, barras deterministicas, caso tabela4; digitalizado, aproximado:
// dados/fig14_lfs.csv) e da sec. 4.6 (A_LFS = 14,9 m2 a 15 h)
const std::map<double, std::pair<double, double>> kFig14 = {{0., {1.684, 14.77}}, {7.5, {1.626, 15.19}},
                                                            {15., {1.512, 14.9}}};

} // namespace

// ------------------------------------------------------------------------------------------------------------------
TEST_CASE("V5: ajuste de van Genuchten a Brooks-Corey (Tabela 1)", "[V5]") {
    const Config cf = Config::CasoPadrao("tabela1");
    const double q = cf.Chuva1();
    for (const Camada *L : {&cf.sup, &cf.inf}) {
        INFO("camada " << L->nome);
        // erro RMS de S_e em 0 < h <= 1,5 m: o do ajuste (scripts/ajuste_vg.py: 0,0229 e 0,0109) com folga de 10%
        double s2 = 0.;
        const int n = 300;
        for (int i = 1; i <= n; ++i) {
            const double h = 1.5 * i / n;
            const double d = L->SeVG(h * cf.gw) - L->SeBC(h);
            s2 += d * d;
        }
        CHECK(std::sqrt(s2 / n) < (L == &cf.sup ? 0.0252 : 0.0120));
        // a condutividade e a de Brooks-Corey (eq. 4): expoente l + 2 + 2/n_BC
        CHECK(L->PotenciaK() == Approx(L == &cf.sup ? 9.2112 : 18.2672).margin(1.e-3));
        // permanente com 600 mm/ano e gradiente unitario: K(S_e) = q da o mesmo S_e (e theta) que Brooks-Corey; o
        // material da biblioteca, com os parametros do problema, tem de dar esse theta (bissecao em p_c)
        const UnsatParams p = ParametrosCamada(cf, *L);
        double lo = 1.e-6, hi = 1.e3;
        for (int it = 0; it < 200; ++it) {
            const double m = std::sqrt(lo * hi);
            (p.hydro.Evaluate(-m, 0.).krw * L->Ks() > q ? lo : hi) = m;
        }
        const double pc = std::sqrt(lo * hi), thetaVG = L->thetaS * p.hydro.Evaluate(-pc, 0.).S;
        const double SeBC = std::pow(q / L->Ks(), 1. / L->PotenciaK());
        CHECK(thetaVG == Approx(L->Theta(SeBC)).margin(1.e-6));
        CHECK(thetaVG == Approx(L == &cf.sup ? 0.2358 : 0.3519).margin(2.e-4)); // valores de ajuste_vg.py
        // a succao desse estado: van Genuchten x Brooks-Corey (1,067 x 1,088 m e 0,991 x 0,918 m): ate 10%
        const double hBC = std::pow(SeBC, -1. / L->nBC) / (L->alphaBC * 100.);
        CHECK(pc / cf.gw == Approx(hBC).epsilon(0.10));
    }
}

TEST_CASE("V5: LFS de uma coluna confinada lateralmente (solucao analitica)", "[V5]") {
    // coluna de 1 m x 10 m com roletes, camada superior nos 2 m de cima, nivel d'agua em z = 5 m (hidrostatico):
    // u_x = 0, sigma_v = integral do peso da mistura, sigma_h = K0 sigma_v com K0 = nu/(1 - nu); o LFS (eq. 8) com
    // sigma' = sigma - S_e p_w tem forma fechada. Verifica o peso com a umidade, o calculo de tensoes, os sinais e chi
    Config cf = Config::CasoPadrao("tabela1");
    cf.nivel = 5., cf.zFace = 5.;
    cf.nu = 0.25; // K0 = 1/3: a parte de baixo da coluna fica com LFS < 1 (a area A_LFS e verificada)
    const double H = 10., zi = 8.;
    auto *g = RectangleMesh(1., H, 1, 40, {}, [&](REAL, REAL z) { return z > zi ? kSuperior : kInferior; });
    Encosta e(cf, g, [](double, double) { return true; });
    e.CondicaoInicial();
    REQUIRE(e.Tensoes());
    // solucao analitica: sigma_v(z) = int_z^H gamma(theta(z')) dz' (regra do ponto medio, 20000 intervalos)
    auto gamma = [&](double z) {
        const Camada &L = z > zi ? cf.sup : cf.inf;
        const double pc = cf.gw * (z - cf.nivel);
        const double theta = pc > 0. ? L.Theta(L.SeVG(pc)) : L.thetaS;
        return ((1. - L.thetaS) * cf.rhoS + theta * 1000.) * cf.g / 1000.;
    };
    auto sigv = [&](double z) {
        const int n = 20000;
        double s = 0.;
        for (int i = 0; i < n; ++i) s += gamma(z + (H - z) * (i + 0.5) / n);
        return s * (H - z) / n;
    };
    const double K0 = cf.nu / (1. - cf.nu), phi = cf.Phi();
    double maxErrS = 0., maxErrK = 0., maxErrL = 0., somaA = 0., somaL = 0., area = 0., areaAn = 0.;
    for (size_t i = 0; i < e.pontos.size(); ++i) {
        const auto &p = e.pontos[i];
        const auto &s = e.estado[i];
        const double sv = sigv(p.z);
        const double pw = cf.gw * (cf.nivel - p.z);
        const Camada &L = p.camada == kSuperior ? cf.sup : cf.inf;
        const double Se = pw >= 0. ? 1. : L.SeVG(-pw);
        const double lfs = std::cos(phi) * (L.c + (0.5 * (1. + K0) * sv - Se * pw) * std::tan(phi)) / (0.5 * (1. - K0) * sv);
        maxErrS = std::max(maxErrS, std::fabs(-s.syy - sv) / (sv + 1.));
        if (sv > 5.) maxErrK = std::max(maxErrK, std::fabs(s.sxx / s.syy - K0));
        if (sv > 5.) maxErrL = std::max(maxErrL, std::fabs(s.lfs / lfs - 1.));
        // sigma_xy = 0 na solucao exata; o P2 com o peso que varia com theta(z) deixa um cisalhamento espurio de
        // ordem 1e-6 sigma_v (1e-4 kPa)
        CHECK(std::fabs(s.sxy) < 1.e-5 * (sv + 10.));
        CHECK(s.Se == Approx(Se).margin(1.e-6));
        area += p.w, somaL += p.w * s.lfs, somaA += p.w * (s.lfs < 1.);
        areaAn += p.w * (lfs < 1.);
    }
    std::cout << "coluna: erro relativo maximo de sigma_v " << maxErrS << ", de K0 " << maxErrK << ", do LFS " << maxErrL
              << "; A(LFS < 1) " << somaA << " m2 (analitico " << areaAn << ")" << std::endl;
    // P2 representa sigma_v quadratico por partes; o peso varia com theta(z): erro de discretizacao pequeno
    // (erros obtidos: 3e-5, 1e-6 e 1e-5; os limites tem folga de ~5 vezes)
    CHECK(maxErrS < 2.e-4);
    CHECK(maxErrK < 1.e-5);
    CHECK(maxErrL < 1.e-4);
    // as metricas da zona (aqui a coluna toda) sao as somas nos pontos de integracao
    const Metricas m = e.CalculaMetricas(0.);
    CHECK(m.areaZona == Approx(10.).epsilon(1.e-9));
    CHECK(m.lfsMedio == Approx(somaL / area).epsilon(1.e-9));
    CHECK(areaAn > 1.); // com nu = 0,25 ha LFS < 1 na parte de baixo
    CHECK(m.aLFS == Approx(areaAn).margin(1.e-9)); // os mesmos pontos de integracao, LFS iguais a 1e-5
}

TEST_CASE("V5: fluxo de 10 anos e chuva de 15 h contra as Figs. 4 e 5 (tabela1, malha grossa)", "[V5]") {
    Rodada &r = Roda("tabela1", true);
    REQUIRE(r.ok1);
    REQUIRE(r.ok2);
    // permanente no fim da etapa 1: a saida pela lateral direita (carga e face) iguala a chuva (balanco < 0,5%)
    const double chuva = r.fim1.chuva / r.fim1.dt, saida = -(r.fim1.dir + r.fim1.face) / r.fim1.dt;
    CHECK(saida == Approx(chuva).epsilon(5.e-3));
    CHECK(std::fabs(r.fim1.residuo / r.fim1.dt) < 1.e-6 * chuva);
    const Instante &i0 = r.inst.at(0.), &i7 = r.inst.at(7.5), &i15 = r.inst.at(15.);
    // Fig. 4, 0 h, crista (x_zona = 2 m), meio da camada superior: theta = 0,236 (digitalizado, 0,236 a 0,239 entre 0,4
    // e 1,6 m de profundidade; o gradiente unitario com Brooks-Corey da 0,2358). Tolerancia 0,01: a cor da Fig. 4 nessa
    // faixa varia 0,01 por pixel de cor
    CHECK(i0.thetaSupCrista == Approx(0.236).margin(0.01));
    // Fig. 5, 0 h, crista, camada superior: tensao de succao 5,4 kPa (digitalizado, 5,35 a 5,58); 0,6 kPa: o ajuste de
    // van Genuchten da h = 1,067 m contra 1,088 m de Brooks-Corey no estado de gradiente unitario (2%), mais a leitura
    CHECK(i0.sigsSupCrista == Approx(5.45).margin(0.6));
    // Fig. 5, 0 h: nivel d'agua (tensao de succao nula) na crista em z_zona ~ 4,1 a 4,4 m (z = 21,1 a 21,4 m).
    // Tolerancia 0,6 m: a carga da lateral direita e o seu nivel nao sao dados (supostos 17 m, como o nivel inicial)
    // e a franja capilar de Brooks-Corey (0,37 m na camada inferior) nao existe no ajuste de van Genuchten
    CHECK(i0.nivelCrista == Approx(21.25).margin(0.6));
    // Fig. 4: frente de umedecimento na crista a ~0,85 m (7,5 h) e ~1,7 m (15 h) (digitalizado: theta passa de 0,32
    // entre 0,83 e 0,9 m e entre 1,6 e 1,8 m). O balanco de massa da 150 mm / (0,402 - 0,236) = 0,90 m e 1,8 m.
    // Tolerancias 0,2 m e 0,25 m: a malha grossa (6 linhas de 0,33 m na camada) espalha a frente
    CHECK(i7.frenteCrista == Approx(0.87).margin(0.2));
    CHECK(i15.frenteCrista == Approx(1.72).margin(0.25));
    // Fig. 4, 15 h: a camada superior esta quase saturada no patamar do pe (theta >= 0,39) e na crista acima da frente
    CHECK(i15.thetaPe > 0.39);
    CHECK(i15.thetaSupCrista > 0.39);
    // mapas inteiros (Figs. 4 e 5, digitalizados) na zona: diferenca RMS de theta <= 0,03 e de tensao de succao
    // <= 2 kPa em cada instante (a frente e a interface, onde o mapa muda de cor em poucos centimetros, dominam a
    // diferenca; a tensao de succao na zona umedecida e ~1 kPa menor com van Genuchten: README.md)
    for (double t : {0., 7.5, 15.}) {
        INFO("t = " << t << " h: theta media " << r.difTheta.at(t).media << " RMS " << r.difTheta.at(t).rms
                    << "; succao media " << r.difSigs.at(t).media << " RMS " << r.difSigs.at(t).rms);
        CHECK(r.difTheta.at(t).n > 1000);
        CHECK(r.difTheta.at(t).rms < 0.03);
        CHECK(r.difSigs.at(t).rms < 2.0);
    }
}

TEST_CASE("V5: LFS medio e area com LFS < 1 contra a Fig. 14 (tabela4, malha grossa)", "[V5]") {
    Rodada &r = Roda("tabela4", true);
    REQUIRE(r.ok1);
    REQUIRE(r.ok2);
    for (const auto &kv : kFig14) {
        const double t = kv.first;
        const Metricas &m = r.inst.at(t).m;
        INFO("t = " << t << " h: LFS medio " << m.lfsMedio << " (limitado a 3: " << m.lfsMedioLim3 << "), A(LFS<1) = "
                    << m.aLFS << " m2; artigo " << kv.second.first << " e " << kv.second.second << " m2");
        CHECK(m.areaZona == Approx(113.).epsilon(1.e-6)); // 4 x 8 + 15 x 5 + 3 x 2 m2 (sec. 4.6: "about 113 m2")
        // media do LFS limitado a 3 (o topo da escala das Figs. 6 e 9): a media literal da eq. (11) inclui pontos com
        // sigma_II -> 0 (LFS -> infinito) logo abaixo da crista e fica ~35% acima da Fig. 14 (README.md). Tolerancia
        // 6%: a diferenca entre a malha grossa e a fina (README.md) mais a leitura das barras (0,005)
        CHECK(m.lfsMedioLim3 == Approx(kv.second.first).epsilon(0.06));
        CHECK(m.lfsMedio > m.lfsMedioLim3);
    }
    // area com LFS < 1 a 0 e 7,5 h: Fig. 14, 14,8 e 15,2 m2. Tolerancia 25%: a area depende do LFS perto de 1 numa
    // regiao larga (a parte funda da zona, sob a crista), sensivel ao nivel d'agua (+-0,3 m contra a Fig. 5) e a malha
    CHECK(r.inst.at(0.).m.aLFS == Approx(kFig14.at(0.).second).epsilon(0.25));
    CHECK(r.inst.at(7.5).m.aLFS == Approx(kFig14.at(7.5).second).epsilon(0.25));
    // 15 h: o artigo da 14,9 m2 (sec. 4.6); aqui a area cresce com a subida do nivel d'agua no pe (ver README.md). Este
    // teste so exige a ordem de grandeza (entre 0,5 e 2 vezes) e o sinal da mudanca de 0 a 15 h (nao diminui)
    CHECK(r.inst.at(15.).m.aLFS > 0.5 * kFig14.at(15.).second);
    CHECK(r.inst.at(15.).m.aLFS < 2.0 * kFig14.at(15.).second);
    CHECK(r.inst.at(15.).m.aLFS >= r.inst.at(0.).m.aLFS - 0.5);
    // o LFS medio cai durante a chuva (Fig. 14: 1,68 -> 1,63 -> 1,51)
    CHECK(r.inst.at(7.5).m.lfsMedioLim3 < r.inst.at(0.).m.lfsMedioLim3);
    CHECK(r.inst.at(15.).m.lfsMedioLim3 < r.inst.at(7.5).m.lfsMedioLim3);
}

TEST_CASE("V5: caso completo na malha fina (tabela4 e tabela1)", "[V5][.lento]") {
    Rodada &r4 = Roda("tabela4", false);
    REQUIRE(r4.ok1);
    REQUIRE(r4.ok2);
    for (const auto &kv : kFig14) {
        const Metricas &m = r4.inst.at(kv.first).m;
        INFO("t = " << kv.first << " h: LFS medio " << m.lfsMedio << " (lim. 3: " << m.lfsMedioLim3 << "), A = "
                    << m.aLFS << " m2");
        CHECK(m.lfsMedioLim3 == Approx(kv.second.first).epsilon(0.05));
    }
    CHECK(r4.inst.at(0.).m.aLFS == Approx(kFig14.at(0.).second).epsilon(0.25));
    CHECK(r4.inst.at(7.5).m.aLFS == Approx(kFig14.at(7.5).second).epsilon(0.25));
    Rodada &r1 = Roda("tabela1", false);
    REQUIRE(r1.ok1);
    REQUIRE(r1.ok2);
    CHECK(r1.inst.at(0.).thetaSupCrista == Approx(0.236).margin(0.01));
    CHECK(r1.inst.at(0.).nivelCrista == Approx(21.25).margin(0.6));
    CHECK(r1.inst.at(7.5).frenteCrista == Approx(0.87).margin(0.15));
    CHECK(r1.inst.at(15.).frenteCrista == Approx(1.72).margin(0.2));
    for (double t : {0., 7.5, 15.}) {
        INFO("t = " << t << " h: theta RMS " << r1.difTheta.at(t).rms << ", succao RMS " << r1.difSigs.at(t).rms
                    << ", LFS (lim. 3) media " << r1.difLFS.at(t).media << " RMS " << r1.difLFS.at(t).rms);
        CHECK(r1.difTheta.at(t).rms < 0.03);
        CHECK(r1.difSigs.at(t).rms < 2.0);
        CHECK(r1.difLFS.at(t).rms < 0.5);
    }
}
