// V6_MediaLocalKL: testes Catch2 (os demorados tem a marca [.lento]).
// O gerador de Karhunen-Loeve de Comum/KLSeparavel.h contra a reducao de variancia da media no arco (Yuan et al. 2021,
// Eqs. 3, 13, 14 e 15) e contra o nucleo exponencial separavel. Tolerancias estatisticas: 4 erros padrao (EP) da
// estatistica amostral (semente fixa: os testes sao deterministicos; com sementes ao acaso, cada comparacao falharia
// com probabilidade ~6e-5).
#define CATCH_CONFIG_MAIN
#include <catch2/catch.hpp>

#include "problema.h"

using namespace v6;

namespace {

/// O caso H = 10 m, beta = 40 graus (Tabela 5) com a KL de 90% no dominio inteiro da Fig. 9 (51,9 x 20 m, M = 11024)
/// e 20000 realizacoes (uns 10 s); calculado uma vez para os testes que o usam
const ResultadoCaso &CasoH10() {
    static const ResultadoCaso R = [] {
        Config c;
        c.caso = "b40_h10";
        c.N = 20000;
        c.nConferir = 20;
        return RodaCaso(c);
    }();
    return R;
}

/// Arcos (Bishop e Fig. 15) e fatores de reducao dos 10 casos (convencao da Eq. 13), calculados uma vez
const std::vector<LinhaTabela> &Tabelas() {
    static const std::vector<LinhaTabela> T = CalculaTabelas(Config{});
    return T;
}

/// int_{-a}^{a} exp(-|s - t| / L) f_i(t) dt por Gauss-Legendre composto, quebrado em t = s (o bico do nucleo)
double IntegralNucleo(const klsep::Modos1D &m, int i, double s) {
    std::vector<double> gx, gw;
    GaussLegendre(10, gx, gw);
    const double h = std::min(m.L, kPi / m.omega[i]) / 4.;
    double tot = 0.;
    for (auto [t0, t1] : {std::pair<double, double>{-m.a, s}, {s, m.a}}) {
        if (t1 - t0 <= 0.) continue;
        const int n = std::max(1, (int)std::ceil((t1 - t0) / h));
        const double dt = (t1 - t0) / n;
        for (int p = 0; p < n; ++p)
            for (int j = 0; j < 10; ++j) {
                const double t = t0 + dt * (p + 0.5 * (gx[j] + 1.));
                tot += 0.5 * dt * gw[j] * std::exp(-std::fabs(s - t) / m.L) * m.Avalia(i, t);
            }
    }
    return tot;
}

/// int_{-a}^{a} f_i f_j dt
double Produto(const klsep::Modos1D &m, int i, int j) {
    std::vector<double> gx, gw;
    GaussLegendre(10, gx, gw);
    const double h = kPi / std::max(m.omega[i], m.omega[j]) / 2.;
    const int n = std::max(1, (int)std::ceil(2. * m.a / h));
    const double dt = 2. * m.a / n;
    double tot = 0.;
    for (int p = 0; p < n; ++p)
        for (int k = 0; k < 10; ++k) {
            const double t = -m.a + dt * (p + 0.5 * (gx[k] + 1.));
            tot += 0.5 * dt * gw[k] * m.Avalia(i, t) * m.Avalia(j, t);
        }
    return tot;
}

} // namespace

TEST_CASE("V6: geometria do arco pelo pe (Eqs. 4 a 12) e os arcos dos casos", "[V6]") {
    // circulo da Fig. 15 (centro digitalizado): as Eqs. 4, 6, 7, 8, 11 e 12 sao identidades geometricas, conferidas com
    // tolerancia de arredondamento (1e-12 relativo)
    const Arco A = Arco::PeloPe(kFig15A, kFig15B, 50.);
    CHECK(A.r == Approx(std::hypot(kFig15A, kFig15B)).epsilon(1.e-14)); // Eq. 4
    CHECK(A.c == Approx(std::sqrt(A.r * A.r - (50. - kFig15B) * (50. - kFig15B)) + kFig15A).epsilon(1.e-14)); // Eq. 6
    // a saida lida na Fig. 15 e (60, 50) m, a meio pixel (0,1 m) de leitura e ~0,5 m de espessura do traco
    CHECK(A.c == Approx(60.).margin(0.5));
    double x, y;
    A.Ponto(0., x, y);
    CHECK(std::fabs(x) < 1.e-12);
    CHECK(std::fabs(y) < 1.e-12);
    A.Ponto(A.T, x, y); // Eq. 12 em l = T = r alpha (Eqs. 7 e 8) cai em A
    CHECK(x == Approx(A.c).epsilon(1.e-12));
    CHECK(y == Approx(50.).epsilon(1.e-12));
    for (int i = 0; i <= 50; ++i) {
        const double l1 = A.T * i / 50., l2 = A.T * ((7 * i) % 51) / 50.;
        double x1, y1, x2, y2;
        A.Ponto(l1, x1, y1);
        A.Ponto(l2, x2, y2);
        CHECK(std::hypot(x1 - A.a, y1 - A.b) == Approx(A.r).epsilon(1.e-12));
        CHECK(std::hypot(x2 - x1, y2 - y1) == Approx(2. * A.r * std::fabs(std::sin((l2 - l1) / (2. * A.r)))).margin(1.e-9));
    }
    CHECK(A.theta == Approx(-std::atan(kFig15A / kFig15B)).epsilon(1.e-14));                 // Eq. 11
    CHECK(std::cos(A.theta) * A.a + std::sin(A.theta) * A.b == Approx(0.).margin(1.e-10)); // tangente em O _|_ OC
    const Talude t{50., 45.};
    CHECK(ArcoDentro(A, t));
    CHECK(A.c > t.XCrista());
    // os 10 arcos (9 circulos de Bishop e o da Fig. 15) ficam no solo e saem na crista
    for (const auto &l : Tabelas()) {
        INFO(l.caso.nome);
        CHECK(ArcoDentro(l.ac.arco, l.caso.t));
        CHECK(l.ac.arco.d == Approx(l.caso.t.H).epsilon(1.e-12));
        CHECK(l.ac.arco.c > l.caso.t.XCrista());
    }
}

TEST_CASE("V6: Eq. 14 contra a forma fechada do segmento reto e a forca bruta; o erro de sinal da Eq. 15", "[V6]") {
    const Correlacao rho = Correlacao::Eq13(kDelta1, kDelta2);
    CHECK(rho.Lx == 1.5); // Eq. 13: exp(-2|tau| / delta) = exp(-|tau| / L) com L = delta / 2
    CHECK(rho.Ly == 0.5);
    // (1) segmento reto (raio de 1e12 m: a curvatura desvia a corda de u / R e entra em |dx| / Lx linearmente, um efeito
    // relativo ~ s T / (R L) < 1e-9; com R = 1e7 m ele ja seria 3e-7): a forma fechada de Vanmarcke com
    // L_e = 1 / (|cos phi| / Lx + |sin phi| / Ly); a quadratura converge exponencialmente: tolerancia 1e-8
    for (double phi : {0., 30., 60., 90.})
        for (double T : {1., 20., 80.}) {
            const double p = phi * kPi / 180., R = 1.e12;
            const Arco A = Arco::ComComprimento(-R * std::sin(p), R * std::cos(p), T);
            const double Le = 1. / (std::fabs(std::cos(p)) / rho.Lx + std::fabs(std::sin(p)) / rho.Ly);
            INFO("phi = " << phi << ", T = " << T);
            CHECK(FatorVariancia(A, rho) == Approx(FatorVarianciaSegmento(T, Le)).epsilon(1.e-8));
        }
    // (2) arcos dos casos: a forca bruta (ponto medio n x n nas coordenadas da Eq. 12, erro O(h^2) pelo bico de
    // exp(-|tau|) na diagonal; Richardson com n e 2n) confere a forma em (s, u) da Eq. 15 corrigida
    const Arco arcos[] = {Arco::PeloPe(kFig15A, kFig15B, 50.), ArcoDoCaso(BuscaCaso("b40_h10")).arco};
    for (const Arco &A : arcos) {
        INFO("arco r = " << A.r << ", T = " << A.T);
        const double g = FatorVariancia(A, rho);
        const double b1 = FatorVarianciaPontoMedio(A, rho, 1500), b2 = FatorVarianciaPontoMedio(A, rho, 3000);
        CHECK(b2 == Approx(g).epsilon(1.e-3));                    // o ponto medio sozinho: erro ~ h^2 / (6 Ly^2)
        CHECK((4. * b2 - b1) / 3. == Approx(g).epsilon(2.e-5));   // extrapolado
        OpcoesQuad fina;
        fina.fator = 0.5;
        fina.ng = 10;
        CHECK(FatorVariancia(A, rho, fina) == Approx(g).epsilon(1.e-10)); // a quadratura ja convergiu
        // (3) a ultima linha da Eq. 15 como impressa (cos((eta + xi)/2r - theta)) NAO e a Eq. 14: difere > 1% no desvio
        const double gi = FatorVarianciaEq15Impressa(A, rho);
        CHECK(std::fabs(std::sqrt(gi / g) - 1.) > 0.01);
    }
}

TEST_CASE("V6: KLSeparavel: autopares 1D do nucleo exponencial (equacao integral e ortonormalidade)", "[V6]") {
    // os autopares 1D (forma fechada, Ghanem e Spanos 1991) sao a base da KL 2D: confere-se a equacao integral
    // int exp(-|s - t| / L) f(t) dt = lambda f(s) e a ortonormalidade por quadratura independente (erro < 1e-10)
    for (auto [D, L] : {std::pair<double, double>{51.92, 1.5}, {20., 0.5}, {259.6, 1.5}}) {
        const auto m = klsep::CalculaModos1D(D, L, 400);
        INFO("D = " << D << ", L = " << L);
        for (int i = 1; i < 400; ++i) CHECK(m.lambda[i] <= m.lambda[i - 1]);
        const int modos[] = {0, 1, 2, 5, 17, 100, 399};
        for (int i : modos) {
            for (double f : {-1., -0.6, -0.1, 0., 0.33, 0.8, 1.}) {
                const double s = f * m.a;
                CHECK(IntegralNucleo(m, i, s) == Approx(m.lambda[i] * m.Avalia(i, s)).margin(1.e-9 * m.lambda[i] * m.norma[i]));
            }
            for (int j : modos) CHECK(Produto(m, i, j) == Approx(i == j ? 1. : 0.).margin(1.e-10));
        }
    }
}

TEST_CASE("V6: KL truncada no talude H = 10 m: fracao retida, variancia no ponto e covariancia", "[V6]") {
    const auto &R = CasoH10();
    // M pelo criterio da Eq. 29 (fracao >= 90%); o passo de um termo e lambda_M / (Dx Dy) ~ 1e-5
    CHECK(R.fracao >= 0.90);
    CHECK(R.fracao < 0.9001);
    CHECK(R.M == 11024);
    // variancia exata do campo truncado (sum lambda_k phi_k^2) nas 20 sondas: a media no retangulo e a fracao retida;
    // os termos que faltam sao de alta frequencia (comprimento de onda < ~0,3 m), entao a variancia local fica a < 1%
    // da fracao longe das bordas (observado: 0,9004 +- 0,0001)
    for (size_t i = 0; i < R.varKL.size(); ++i) {
        INFO(R.sondas.nome[i]);
        CHECK(R.varKL[i] == Approx(R.fracao).margin(0.01));
    }
    // covariancia exata do campo truncado contra o nucleo: |C_M - C| <= sqrt((1 - v_P)(1 - v_Q)) (Cauchy-Schwarz nos
    // termos que faltam: um limite rigoroso, ~0,1 aqui)
    for (size_t p = 0; p < R.sondas.pares.size(); ++p) {
        const int i = R.sondas.pares[p].first, j = R.sondas.pares[p].second;
        INFO(R.sondas.nome[j]);
        const double CM = R.rhoKL[p] * std::sqrt(R.varKL[i] * R.varKL[j]);
        CHECK(std::fabs(CM - R.rhoNucleo[p]) <= std::sqrt((1. - R.varKL[i]) * (1. - R.varKL[j])));
        // a parte que falta e de escala curta: em defasagens >= 1 L (x e y: pares 2 a 4 e 7 a 9) a covariancia
        // truncada fica a menos de 3% do nucleo (observado: -2,0%), contra 10% em defasagem 0 (a variancia)
        if ((p >= 2 && p <= 4) || (p >= 7 && p <= 9)) CHECK(CM == Approx(R.rhoNucleo[p]).epsilon(0.03));
    }
    // a quadratura do arco usada para as medias dos termos F_k ja convergiu: com paineis duas vezes menores a variancia
    // exata da media no arco muda menos de 1e-8 (relativo)
    const Talude &t = R.caso.t;
    const klsep::KL2D kl(t.X0(), t.Y0(), t.X1() - t.X0(), t.Y1() - t.Y0(), R.rho.Lx, R.rho.Ly, 0.90, 400000);
    REQUIRE(kl.M() == R.M);
    const QuadArco q = CriaQuadArcoParaKL(R.ac.arco, kl);
    CHECK(q.N() == R.nq);
    const QuadArco q2 = CriaQuadArco(R.ac.arco, R.ac.arco.T / (q.N() / 8) / 2., 8);
    const double V1 = VarianciaMediaAcumulada(kl, MediasDosTermos(kl, q)).back();
    const double V2 = VarianciaMediaAcumulada(kl, MediasDosTermos(kl, q2)).back();
    CHECK(V1 == Approx(R.G2KL).epsilon(1.e-14));
    CHECK(V2 == Approx(V1).epsilon(1.e-8));
}

TEST_CASE("V6: Monte Carlo no talude H = 10 m: variancia no ponto, correlacao nas defasagens e media no arco", "[V6]") {
    const auto &R = CasoH10();
    const auto &mc = R.mc;
    const long N = mc.N;
    REQUIRE(N == 20000);
    // (a) a media no arco por KL2D::Campo nos 1480 pontos de quadratura e por sum sqrt(lambda_k) xi_k F_k e a mesma
    // conta reordenada: difere so por arredondamento
    CHECK(mc.nConferidas == 20);
    CHECK(mc.maxDifCampo < 1.e-12);
    // (b) variancia e media amostrais nas sondas contra a variancia exata do campo truncado (4 EP)
    for (size_t i = 0; i < R.sondas.x.size(); ++i) {
        INFO(R.sondas.nome[i]);
        CHECK(std::fabs(mc.sonda[i].Var() - R.varKL[i]) < 4. * EPdaVar(R.varKL[i], N));
        CHECK(std::fabs(mc.sonda[i].media) < 4. * std::sqrt(R.varKL[i] / N));
    }
    // (c) correlacao amostral nas defasagens 0,25 a 3 L (x, y, diagonais e ao longo do arco) contra a correlacao exata
    // do campo truncado (4 EP) e contra o nucleo (o limite de Cauchy-Schwarz + 4 EP)
    for (size_t p = 0; p < R.sondas.pares.size(); ++p) {
        const int i = R.sondas.pares[p].first, j = R.sondas.pares[p].second;
        INFO(R.sondas.nome[j] << ": nucleo " << R.rhoNucleo[p] << ", KL " << R.rhoKL[p] << ", MC " << mc.par[p].Rho());
        const double ep = EPdaCorrelacao(R.rhoKL[p], N);
        CHECK(std::fabs(mc.par[p].Rho() - R.rhoKL[p]) < 4. * ep);
        CHECK(std::fabs(mc.par[p].Rho() * std::sqrt(R.varKL[i] * R.varKL[j]) - R.rhoNucleo[p]) <
              std::sqrt((1. - R.varKL[i]) * (1. - R.varKL[j])) + 4. * ep);
    }
    // (d) media no arco: media 0, desvio contra o exato da KL truncada (4 EP: EP = dp / sqrt(2 (N - 1)) = 0,5%) e contra a
    // Eq. 14 (a diferenca exata Gamma - Gamma_KL, 0,3%, mais 4 EP); normalidade (assimetria, curtose e cauda de 5%)
    const double G = std::sqrt(R.G2), GKL = std::sqrt(R.G2KL), sd = mc.media.DP(), ep = EPdoDP(sd, N);
    INFO("Gamma = " << G << ", Gamma_KL = " << GKL << ", dp MC = " << sd << " +- " << ep);
    CHECK(mc.dpExato == Approx(GKL).epsilon(1.e-12));
    CHECK(std::fabs(mc.media.media) < 4. * GKL / std::sqrt((double)N));
    CHECK(std::fabs(sd - GKL) < 4. * ep);
    CHECK(std::fabs(sd - G) < (G - GKL) + 4. * ep);
    CHECK(std::fabs(mc.Assimetria()) < 4. * std::sqrt(15. / N));
    CHECK(std::fabs(mc.ExcessoCurtose()) < 4. * std::sqrt(96. / N));
    CHECK(std::fabs((double)mc.nFora / N - 0.05) < 4. * std::sqrt(0.05 * 0.95 / N));
}

TEST_CASE("V6: convergencia com o numero de termos da KL (talude H = 10 m)", "[V6]") {
    const auto &R = CasoH10();
    const double G = std::sqrt(R.G2);
    REQUIRE(R.conv.size() >= 6);
    for (size_t i = 0; i < R.conv.size(); ++i) {
        const auto &c = R.conv[i];
        INFO("M = " << c.M << ": fracao " << c.fracao << ", Gamma_KL / Gamma " << std::sqrt(c.G2) / G);
        // Mercer: Var[X_M] = sum_{k<M} lambda_k F_k^2 cresce com M e nunca passa de Gamma^2 (erro de quadratura ~1e-7)
        CHECK(c.G2 <= R.G2 * (1. + 1.e-6));
        CHECK(c.varP0 <= 1. + 1.e-9);
        if (i > 0) {
            CHECK(c.G2 >= R.conv[i - 1].G2);
            CHECK(c.fracao >= R.conv[i - 1].fracao);
            CHECK(c.varP0 >= R.conv[i - 1].varP0 - 1.e-3); // a variancia num ponto cresce na media (pode oscilar pouco)
        }
        // a correlacao normalizada do campo truncado em 1 L fica acima do nucleo (e^-1) e desce para ele com M
        CHECK(c.rhoLx > std::exp(-1.));
        CHECK(c.rhoLy > std::exp(-1.));
    }
    // a media no arco converge muito mais depressa que a variancia no ponto: com M = 1000 (62% da variancia) o desvio
    // ja esta a 4,1% de Gamma; com 90% (Eq. 29), a 0,3%. O deficit de variancia da media (1 - Gamma_KL^2 / Gamma^2) e
    // < 15% do deficit da variancia no ponto (1 - fracao): observado 6%
    for (const auto &c : R.conv)
        if (c.M == 1000) CHECK(std::sqrt(c.G2) / G > 0.95);
    CHECK(std::sqrt(R.G2KL) / G > 0.99);
    CHECK(1. - R.G2KL / R.G2 < 0.15 * (1. - R.fracao));
    // e a correlacao normalizada do campo truncado em 1 L e C_M / v_M ~ e^-1 / fracao (a covariancia em 1 L fica a < 3%
    // do nucleo, a variancia cai 10%): 0,400 contra 0,368 do nucleo (o campo truncado e mais suave)
    CHECK(R.conv.back().rhoLx == Approx(std::exp(-1.) / R.fracao).epsilon(0.03));
    CHECK(R.conv.back().rhoLy == Approx(std::exp(-1.) / R.fracao).epsilon(0.03));
}

TEST_CASE("V6: arcos reconstruidos contra as Tabelas 5 a 8 (Bishop com gamma calibrado e o circulo da Fig. 15)", "[V6]") {
    const auto &T = Tabelas();
    // gamma calibrado: FS de Bishop = 0,9938 no talude H = 50 m, beta = 40 graus (Tabela 7); valor plausivel para loess
    const double g = GammaCalibrado();
    CHECK(g > 14.);
    CHECK(g < 22.);
    for (const auto &l : T) {
        INFO(l.caso.nome << ": FS Bishop " << l.ac.FSbishop << " (artigo " << l.caso.FS << "), Gamma " << l.G
                         << " (artigo c " << l.caso.scL / kSigmaC << ", phi " << l.caso.sphiL / kSigmaPhi << ")");
        if (l.caso.nome == "b40_h50") CHECK(l.ac.FSbishop == Approx(kFSCalibracao).epsilon(1.e-6));
        // os outros FS sao previsoes com o gamma calibrado: Bishop simplificado contra o SRM do ABAQUS com psi = 0.
        // Tolerancia de 10%: o proprio artigo difere ate ~10% entre os FS impressos e os das curvas de estado limite
        // ajustadas (Tabelas 3, 4, 7 e 8); o circulo da Fig. 15 nao e critico para Bishop (FS maior que o minimo)
        if (!l.caso.fig15) CHECK(l.ac.FSbishop == Approx(l.caso.FS).epsilon(0.10));
        // Gamma contra as Tabelas 5 e 6. NAO ha concordancia: o nosso Gamma (Eq. 14 integrada) fica de 0 a 20% abaixo
        // das tabelas em 18 dos 20 valores (+3% no maximo). A faixa abaixo apenas delimita a discrepancia documentada
        // no README (detecta, por exemplo, uma troca de convencao de delta, que da +40%):
        //  - os arcos do artigo nao sao impressos (os nossos sao reconstruidos);
        //  - a reducao de c e de phi tem de ser a mesma para o mesmo arco (o fator nao depende de sigma), mas as tabelas
        //    dao fatores que diferem ate 9% (precisao das proprias tabelas);
        //  - com o circulo do proprio artigo (Fig. 15) o nosso valor e 19% (c) e 13% (phi) menor que a Tabela 6.
        const double rc = l.G / (l.caso.scL / kSigmaC), rp = l.G / (l.caso.sphiL / kSigmaPhi);
        CHECK(rc > 0.78);
        CHECK(rc < 1.06);
        CHECK(rp > 0.78);
        CHECK(rp < 1.06);
        // a leitura L = delta (secao 4 do artigo) da sempre um Gamma ~40% maior
        CHECK(l.Goutra / l.G > 1.3);
    }
    // a inconsistencia interna das tabelas (fatores de c e de phi para o mesmo arco): ate 9%
    double maxdif = 0.;
    for (const auto &l : T) maxdif = std::max(maxdif, std::fabs((l.caso.scL / kSigmaC) / (l.caso.sphiL / kSigmaPhi) - 1.));
    CHECK(maxdif > 0.08);
    // Fig. 15 (H = 50 m, beta = 45 graus) contra o nosso circulo critico de Bishop do mesmo talude: raios diferentes
    // (125 contra 83 m) mas comprimentos e fatores proximos: o fator depende sobretudo do comprimento do arco
    const LinhaTabela *f15 = nullptr, *b45 = nullptr;
    for (const auto &l : T) {
        if (l.caso.nome == "fig15") f15 = &l;
        if (l.caso.nome == "h50_b45") b45 = &l;
    }
    REQUIRE(f15);
    REQUIRE(b45);
    CHECK(b45->ac.arco.T == Approx(f15->ac.arco.T).epsilon(0.05));
    CHECK(b45->G == Approx(f15->G).epsilon(0.02));
}

// ------------------------------------------------------------------------------------------------ [.lento]

namespace {
void ConfereCasoGrande(const ResultadoCaso &R, double minGKL) {
    const auto &mc = R.mc;
    const long N = mc.N;
    const double G = std::sqrt(R.G2), GKL = std::sqrt(R.G2KL), sd = mc.media.DP(), ep = EPdoDP(sd, N);
    INFO(R.caso.nome << ": M = " << R.M << ", fracao " << R.fracao << ", Gamma = " << G << ", Gamma_KL = " << GKL
                     << ", dp MC = " << sd << " +- " << ep);
    CHECK(GKL / G > minGKL);
    CHECK(GKL <= G * (1. + 1.e-6));
    CHECK(mc.maxDifCampo < 1.e-12);
    CHECK(std::fabs(sd - GKL) < 4. * ep);
    CHECK(std::fabs(sd - G) < (G - GKL) + 4. * ep);
    CHECK(std::fabs(mc.Assimetria()) < 4. * std::sqrt(15. / N));
    CHECK(std::fabs(mc.ExcessoCurtose()) < 4. * std::sqrt(96. / N));
    for (size_t i = 0; i < R.sondas.x.size(); ++i) {
        INFO(R.sondas.nome[i]);
        CHECK(R.varKL[i] == Approx(R.fracao).margin(0.01));
        CHECK(std::fabs(mc.sonda[i].Var() - R.varKL[i]) < 4. * EPdaVar(R.varKL[i], N));
    }
    for (size_t p = 0; p < R.sondas.pares.size(); ++p) {
        INFO(R.sondas.nome[R.sondas.pares[p].second]);
        CHECK(std::fabs(mc.par[p].Rho() - R.rhoKL[p]) < 4. * EPdaCorrelacao(R.rhoKL[p], N));
    }
}
} // namespace

TEST_CASE("V6: H = 50 m, beta = 40 graus: KL de 90% no dominio inteiro (M = 272700) e 10^4 realizacoes", "[V6][.lento]") {
    Config c;
    c.caso = "b40_h50";
    c.N = 10000;
    c.nConferir = 3;
    const auto R = RodaCaso(c);
    CHECK(R.M == 272700);
    ConfereCasoGrande(R, 0.99);
}

TEST_CASE("V6: arco da Fig. 15: KL de 90% (10^4 realizacoes) e de 95% (sem Monte Carlo)", "[V6][.lento]") {
    Config c;
    c.caso = "fig15";
    c.N = 10000;
    c.nConferir = 3;
    const auto R = RodaCaso(c);
    ConfereCasoGrande(R, 0.99);
    // com 95% da variancia (M ~ 7,5e5) a media no arco fica ainda mais perto da Eq. 14
    Config c95 = c;
    c95.N = 0;
    c95.fracao = 0.95;
    c95.Mmax = 900000;
    const auto R95 = RodaCaso(c95);
    INFO("M95 = " << R95.M << ", Gamma_KL / Gamma = " << std::sqrt(R95.G2KL / R95.G2));
    CHECK(R95.fracao >= 0.95);
    CHECK(R95.G2KL > R.G2KL);
    CHECK(R95.G2KL <= R95.G2 * (1. + 1.e-6));
    CHECK(1. - R95.G2KL / R95.G2 < 0.6 * (1. - R.G2KL / R.G2));
}

TEST_CASE("V6: a outra leitura de delta (L = delta, secao 4 do artigo) no arco da Fig. 15", "[V6][.lento]") {
    Config c;
    c.caso = "fig15";
    c.convencao = "delta";
    c.N = 10000;
    c.nConferir = 3;
    const auto R = RodaCaso(c);
    CHECK(R.rho.Lx == 3.);
    CHECK(R.rho.Ly == 1.);
    ConfereCasoGrande(R, 0.99);
    // Gamma ~ 0,169 contra 0,120 na leitura da Eq. 13: sigma_cL = 1,47 kPa contra 1,28 kPa da Tabela 6
    CHECK(std::sqrt(R.G2) == Approx(0.1688).epsilon(0.01));
}
