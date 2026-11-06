// V1_LiakopoulosTrifasico: testes Catch2. Os padrao rodam os casos ate 30 min (cerca de 45 s cada, 3 casos); os de
// 600 min tem a marca [.lento]. Os valores de referencia vem de dados/*.csv (curvas digitalizadas, aproximadas) e da
// tabela do LearningLiakopoulos (branch UnsaturatedSoil, validada contra a referencia 1D e Callari e Abati 2009).
#define CATCH_CONFIG_MAIN
#include <catch2/catch.hpp>

#include "problema.h"

#include <map>

using namespace v1;

namespace {

/// Rodadas guardadas entre os TEST_CASEs (cada caso roda uma vez por tempo final)
const Resultado &Rodada(const std::string &caso, REAL tFim)
{
    static std::map<std::pair<std::string, REAL>, Resultado> cache;
    auto chave = std::make_pair(caso, tFim);
    auto it = cache.find(chave);
    if (it != cache.end()) return it->second;
    Opcoes op;
    op.tFim = tFim;
    op.verboso = false;
    Resultado r = RodaCaso(CasoPorNome(caso), op);
    std::cout << "[V1] " << caso << " ate " << tFim << " s: " << r.passos << " passos, " << r.seg << " s"
              << std::endl;
    return cache.emplace(chave, std::move(r)).first->second;
}

/// Perfil gravado no instante t (falha o teste se nao existir)
const Perfil &PerfilEm(const Resultado &r, REAL t)
{
    const Perfil *p = r.PerfilEm(t);
    REQUIRE(p != nullptr);
    return *p;
}

REAL Interp(const std::vector<REAL> &x, const std::vector<REAL> &y, REAL xi)
{
    if (xi <= x.front()) return y.front();
    for (size_t i = 1; i < x.size(); ++i)
        if (xi <= x[i]) return y[i - 1] + (y[i] - y[i - 1]) * (xi - x[i - 1]) / (x[i] - x[i - 1]);
    return y.back();
}

} // namespace

// ---------------------------------------------------------------------------------------------------------------
// Leis do ar: o mapeamento dos dados de Laloui et al. nas leis da biblioteca
// ---------------------------------------------------------------------------------------------------------------
TEST_CASE("V1: k_ra = (1 - S_e)^3 com S_r^r = 1% e k_ra >= 1% (Eq. 61 e Secao 5.2 de Laloui et al.)", "[V1]")
{
    // A Eq. (61) e Brooks-Corey (1 - S_e)^2 (1 - S_e^((2 + lambda)/lambda)) com lambda -> infinito; com lambda = 1e6 o
    // erro relativo e da ordem de 2/lambda |ln S_e| < 1e-5. Tolerancia 1e-5 (relativa).
    const UnsatParams p = ParamsTrifasico();
    const TPZUnsaturatedHydraulics &h = p.hydro;
    for (REAL S : {0.05, 0.2, 0.5, 0.7, 0.75}) {
        const REAL Se = (S - 0.01) / 0.99, esperado = std::pow(1. - Se, 3);
        REAL kr, dkr;
        h.AirPermeability(S, kr, dkr);
        INFO("S_w = " << S << ": k_ra = " << kr << ", (1 - S_e)^3 = " << esperado);
        CHECK(kr == Approx(esperado).epsilon(1.e-5));
        CHECK(dkr == Approx(-3. * std::pow(1. - Se, 2) / 0.99).epsilon(1.e-4));
    }
    // acima de S_w = 0,7945 ((1 - S_e)^3 = 0,01) vale o piso de 1%: na coluna (S_w > 0,9) k_ra = 0,01 em toda parte
    for (REAL S : {0.8, 0.9, 0.95, 1.}) {
        REAL kr, dkr;
        h.AirPermeability(S, kr, dkr);
        CHECK(kr == Approx(0.01));
        CHECK(dkr == 0.);
    }
    // mobilidade do ar / mobilidade da agua = (mu_w / mu_a) k_ra / k_rw = 55,6 x 0,01 = 0,556 com a agua saturada
    CHECK(p.muRatio * 0.01 == Approx(0.5556).epsilon(1.e-3));
}

// ---------------------------------------------------------------------------------------------------------------
// Ar passivo: o porte do LearningLiakopoulos (regressao da validacao)
// ---------------------------------------------------------------------------------------------------------------
TEST_CASE("V1: ar passivo reproduz o LearningLiakopoulos validado (referencia 1D e Callari e Abati 2009)", "[V1]")
{
    const Resultado &r = Rodada("bifasico", 1800.);
    REQUIRE(r.ok);
    // Tabela do README do LearningLiakopoulos (NeoPZ com os passos da referencia 1D, igual a referencia 1D nos digitos
    // impressos): t (s), q_b (1e-6 m/s), p_w em y = 0,5 m e no topo (kPa), w do topo (mm). Mesma malha, mesmos passos:
    // a diferenca tolerada e meia unidade do ultimo digito impresso (0,005 kPa, 0,0005 mm, 0,0005e-6 m/s) mais 0,001e-6
    // m/s na vazao (a referencia 1D difere do NeoPZ em 0,001e-6 m/s em 600 s).
    struct Ref {
        REAL t, qb, pwMeio, pwTopo, w;
    };
    const std::vector<Ref> tab = {{300., 3.373, -1.30, -4.15, -0.550},
                                  {600., 2.727, -1.98, -5.22, -0.774},
                                  {1200., 2.011, -2.74, -6.41, -1.020},
                                  {1800., 1.578, -3.20, -7.13, -1.167}};
    for (const Ref &f : tab) {
        const REAL qb = r.HistEm(&Linha::qb, f.t) * 1.e6, pm = r.HistEm(&Linha::pwMeio, f.t);
        const REAL pt = r.HistEm(&Linha::pwTopo, f.t), w = r.HistEm(&Linha::wTopo, f.t) * 1.e3;
        INFO("t = " << f.t << " s: q_b = " << qb << " [" << f.qb << "], p_w(0,5) = " << pm << " [" << f.pwMeio
                    << "], p_w(topo) = " << pt << " [" << f.pwTopo << "], w = " << w << " [" << f.w << "]");
        CHECK(std::fabs(qb - f.qb) <= 0.0015);
        CHECK(std::fabs(pm - f.pwMeio) <= 0.0051);
        CHECK(std::fabs(pt - f.pwTopo) <= 0.0051);
        CHECK(std::fabs(w - f.w) <= 0.00051);
    }
    CHECK(r.pontosPlasticos == 0); // esqueleto elastico (c' = 1e5 kPa)

    // Contra as curvas digitalizadas de Callari e Abati (2009) (dados/callari2009_fig3_fig4_digitalizado.csv,
    // incerteza de leitura ~1% da faixa de cada eixo): o LearningLiakopoulos obteve media de p_w - p_w,artigo de
    // -0,061 kPa nos perfis, media |dq_b|/q_b de 2,0% para t < 2400 s e media |dw| de 0,010 mm. Tolerancias: 0,1 kPa
    // (1% da faixa de 0 a -10 kPa mais a diferenca obtida), 3% e 0,03 mm.
    const std::vector<Curva> cs = LeCurvas(Dados("callari2009_fig3_fig4_digitalizado.csv"));
    REQUIRE(!cs.empty());
    REAL soma = 0.;
    int n = 0;
    for (REAL t : {600., 1200., 1800.}) {
        const Curva *c = AchaCurva(cs, "3b", "p_w em t=" + std::to_string(int(t)) + " s");
        REQUIRE(c != nullptr);
        const Perfil &p = PerfilEm(r, t);
        for (size_t i = 0; i < c->x.size(); ++i, ++n) soma += p.Em(p.pw, c->y[i]) - c->x[i];
    }
    const REAL mediaPw = soma / n;
    soma = 0., n = 0;
    REAL somaW = 0.;
    int nW = 0;
    for (const Curva &c : cs) {
        if (c.figura == "3a")
            for (size_t i = 0; i < c.x.size(); ++i)
                if (c.x[i] > 0. && c.x[i] <= 1800.) {
                    soma += std::fabs(r.HistEm(&Linha::qb, c.x[i]) * 1.e6 - c.y[i]) / c.y[i];
                    n++;
                }
        if (c.figura == "4a" && c.curva == "recalque do topo (parcialmente saturada)")
            for (size_t i = 0; i < c.x.size(); ++i)
                if (c.x[i] <= 1800.) {
                    somaW += std::fabs(r.HistEm(&Linha::wTopo, c.x[i]) * 1.e3 - c.y[i]);
                    nW++;
                }
    }
    INFO("Callari e Abati: media dp_w = " << mediaPw << " kPa, media |dq_b|/q_b = " << soma / n
                                          << ", media |dw| = " << somaW / nW << " mm");
    CHECK(std::fabs(mediaPw) <= 0.1);
    CHECK(soma / n <= 0.03);
    CHECK(somaW / nW <= 0.03);
}

// ---------------------------------------------------------------------------------------------------------------
// Ar ativo com a lei do ar do benchmark: a faixa das cinco equipes (Fig. 6a) e a diferenca duas x tres fases (Fig. 8b)
// ---------------------------------------------------------------------------------------------------------------
TEST_CASE("V1: ar ativo (lei do ar do benchmark): p_a a 30 min dentro da faixa das cinco equipes (Fig. 6a)", "[V1]")
{
    const Resultado &r = Rodada("trifasico_alert", 1800.);
    REQUIRE(r.ok);
    const Perfil &p = PerfilEm(r, 1800.);
    const std::vector<Curva> cs = LeCurvas(Dados("laloui2003_fig6a_pa_30min.csv"));
    const Curva *emin = AchaCurva(cs, "6a", "envoltoria min"), *emax = AchaCurva(cs, "6a", "envoltoria max");
    REQUIRE(emin != nullptr);
    REQUIRE(emax != nullptr);
    // perfil inteiro dentro da envoltoria das cinco curvas, com folga de 0,15 kPa (largura dos marcadores e leitura)
    for (size_t i = 0; i < emin->y.size(); ++i) {
        const REAL y = emin->y[i];
        if (y < 0.04 || y > 0.96) continue; // nos extremos p_a = 0 (Dirichlet) e a envoltoria e so a leitura
        const REAL lo = emin->x[i], hi = Interp(emax->y, emax->x, y), pa = p.Em(p.pa, y);
        INFO("y = " << y << " m: p_a = " << pa << " kPa, faixa [" << lo << ", " << hi << "]");
        CHECK(pa >= lo - 0.15);
        CHECK(pa <= hi + 0.15);
    }
    // minimo: entre o menor e o maior dos minimos das equipes (-6,95 kPa LCPC a -5,28 kPa UPD, folga 0,15 kPa) e na
    // faixa de alturas dos minimos (0,80 a 0,90 m, folga de um elemento, 0,025 m)
    REAL mmin = 0., mmax = -1.e9, ymin = 1., ymax = 0.;
    for (const Curva &c : cs)
        if (c.curva.rfind("minimo", 0) == 0) {
            mmin = std::min(mmin, c.x[0]), mmax = std::max(mmax, c.x[0]);
            ymin = std::min(ymin, c.y[0]), ymax = std::max(ymax, c.y[0]);
        }
    const auto [pa, y] = p.Minimo(p.pa);
    INFO("minimo de p_a em 30 min: " << pa << " kPa em y = " << y << " m; equipes: [" << mmin << ", " << mmax
                                     << "] kPa, [" << ymin << ", " << ymax << "] m");
    CHECK(pa >= mmin - 0.15);
    CHECK(pa <= mmax + 0.15);
    CHECK(y >= ymin - 0.025);
    CHECK(y <= ymax + 0.025);
}

TEST_CASE("V1: ar ativo (lei do ar do benchmark): p_w no topo e diferenca duas x tres fases (Figs. 5a e 8b)", "[V1]")
{
    const Resultado &r2 = Rodada("bifasico", 1800.);
    const Resultado &r3 = Rodada("trifasico_alert", 1800.);
    REQUIRE(r2.ok);
    REQUIRE(r3.ok);
    const std::vector<Curva> cs = LeCurvas(Dados("laloui2003_fig8b_pw_topo.csv"));
    const Curva *c2 = AchaCurva(cs, "8b", "duas fases"), *c3 = AchaCurva(cs, "8b", "tres fases");
    REQUIRE(c2 != nullptr);
    REQUIRE(c3 != nullptr);
    for (REAL tmin : {5., 10., 30.}) {
        const REAL t = 60. * tmin;
        const REAL p2 = PerfilEm(r2, t).pw.back(), p3 = PerfilEm(r3, t).pw.back();
        const REAL a2 = Interp(c2->x, c2->y, tmin), a3 = Interp(c3->x, c3->y, tmin);
        INFO("t = " << tmin << " min: p_w(topo) duas fases " << p2 << " [" << a2 << "], tres fases " << p3 << " ["
                    << a3 << "], diferenca " << p3 - p2 << " [" << a3 - a2 << "]");
        // o ar ativo drena a agua mais depressa no inicio: p_w mais negativo no topo
        CHECK(p3 < p2);
        // diferenca tres - duas fases: tolerancia 0,8 kPa. O modelo de duas fases de Laloui usa a retencao de Seker
        // com S_r <= 99% e o nosso a de Liakopoulos: o nosso p_w de duas fases no topo fica ate 0,63 kPa mais negativo
        // que o dele (30 min); mais 0,1 kPa de leitura de cada curva.
        CHECK(std::fabs((p3 - p2) - (a3 - a2)) <= 0.8);
        // p_w de tres fases no topo contra o do artigo: tolerancia 0,5 kPa (as Figs. 5a e 8b do proprio artigo diferem
        // em ate 0,4 kPa em 30 min: -8,96 e -8,56 kPa)
        CHECK(std::fabs(p3 - a3) <= 0.5);
    }
}

// ---------------------------------------------------------------------------------------------------------------
// Ar ativo com os dados impressos de Laloui et al.: discrepancia documentada
// ---------------------------------------------------------------------------------------------------------------
TEST_CASE("V1: ar ativo com os dados impressos (k_ra >= 1%): p_a a 30 min fica ~10x abaixo da faixa da Fig. 6a "
          "(falha conhecida)",
          "[V1][!shouldfail]")
{
    // Com os dados de Laloui et al. como impressos (k_ra = (1 - S_e)^3 >= 1%, mu_a = 1,8e-5 Pa s, lei de Darcy do ar
    // com a permeabilidade intrinseca k: mobilidade do ar = 0,556 da da agua saturada), o ar entra facilmente pelas
    // duas extremidades e o minimo de p_a em 30 min e cerca de -0,6 kPa, contra -5,3 a -7,0 kPa das cinco equipes.
    // Este teste afirma a concordancia com a Fig. 6a e falha: [!shouldfail] o mantem no ctest como discrepancia
    // registrada (se passar a concordar, o teste acusa). Ver o README (discussao da Fig. 6a).
    const Resultado &r = Rodada("trifasico", 1800.);
    REQUIRE(r.ok);
    const Perfil &p = PerfilEm(r, 1800.);
    const auto [pa, y] = p.Minimo(p.pa);
    const Resultado &r2 = Rodada("bifasico", 1800.);
    const REAL dpw = p.pw.back() - PerfilEm(r2, 1800.).pw.back();
    WARN("trifasico (dados impressos): minimo de p_a em 30 min = " << pa << " kPa em y = " << y
                                                                   << " m (faixa das equipes: -6,95 a -5,28 kPa); "
                                                                      "p_w(topo) tres - duas fases = "
                                                                   << dpw << " kPa (Fig. 8b: -2,1 kPa)");
    CHECK(pa >= -7.1);
    CHECK(pa <= -5.13);
    CHECK(std::fabs(dpw - (-2.06)) <= 0.8);
}

// ---------------------------------------------------------------------------------------------------------------
// Casos completos (600 min): [.lento]
// ---------------------------------------------------------------------------------------------------------------
TEST_CASE("V1: 600 min: minimos de p_a no tempo (Fig. 9a) e p_w hidrostatico no fim (Fig. 5a)", "[V1][.lento]")
{
    const Resultado &r3 = Rodada("trifasico_alert", 36000.);
    REQUIRE(r3.ok);
    const std::vector<Curva> cs = LeCurvas(Dados("laloui2003_fig9a_pa.csv"));
    for (REAL tmin : {5., 10., 30., 60., 120., 600.}) {
        const Curva *c = AchaCurva(cs, "9a", "minimo t=" + std::to_string(int(tmin)) + " min");
        REQUIRE(c != nullptr);
        const Perfil &p = PerfilEm(r3, 60. * tmin);
        const auto [pa, y] = p.Minimo(p.pa);
        INFO("t = " << tmin << " min: minimo de p_a " << pa << " kPa em y = " << y << " m [Fig. 9a: " << c->x[0]
                    << " kPa em " << c->y[0] << " m]");
        // tolerancia 1,7 kPa: a dispersao dos minimos das cinco equipes em 30 min (Fig. 6a, -6,95 a -5,28 kPa); as
        // proprias Figs. 6a e 9a do EPFL diferem em 0,76 kPa em 30 min (-6,04 e -6,80 kPa). Altura: 0,07 m (metade
        // da dispersao das alturas dos minimos, 0,10 m, mais a leitura, 0,02 m)
        CHECK(std::fabs(pa - c->x[0]) <= 1.7);
        CHECK(std::fabs(y - c->y[0]) <= 0.07);
    }
    // 600 min: p_w quase hidrostatico (experimento da Fig. 5a e os dois modelos); tolerancia 0,4 kPa contra o
    // experimento (dispersao dos marcadores em torno da hidrostatica: 0,33 kPa em y = 0,8 m) e 0,1 kPa contra -gamma_w h
    const std::vector<Curva> c5 = LeCurvas(Dados("laloui2003_fig5.csv"));
    const Curva *e600 = AchaCurva(c5, "5a", "experimento t=600 min");
    REQUIRE(e600 != nullptr);
    const Resultado &r2 = Rodada("bifasico", 36000.);
    REQUIRE(r2.ok);
    for (const Resultado *r : {&r2, &r3}) {
        const Perfil &p = PerfilEm(*r, 36000.);
        for (size_t i = 0; i < e600->x.size(); ++i) {
            const REAL y = e600->x[i], pw = p.Em(p.pw, y);
            INFO(r->caso << " 600 min, y = " << y << ": p_w = " << pw << " [experimento " << e600->y[i] << "]");
            CHECK(std::fabs(pw - e600->y[i]) <= 0.4);
            CHECK(std::fabs(pw + 9.81 * y) <= 0.1);
        }
    }
}
