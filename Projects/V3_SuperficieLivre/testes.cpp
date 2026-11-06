// V3_SuperficieLivre: testes Catch2 contra Sharma, Fujisawa & Murakami (2021), FEAD 197:103632 (os demorados tem a
// marca [.lento]: rode "V3_testes [.lento]"). Os testes padrao usam malhas reduzidas (1 m na barragem retangular,
// 0,25 m na vertical da barragem inclinada) e cada caso roda uma so vez por execucao (os resultados ficam guardados).
#define CATCH_CONFIG_MAIN
#include <catch2/catch.hpp>

#include "problema.h"

using namespace v3;

namespace {

/// parametros dos testes padrao: as malhas reduzidas e o regime permanente com TOL_EST = 1e-4
Numerica NumTeste(const Caso &c) {
    Numerica num;
    if (c.nome == "retangular") num.nx = 16, num.ny = 24;   // 1 m x 1 m
    if (c.nome == "inclinada") num.nx = 14, num.ny = 20;    // 0,25 m na vertical, 0,5 m na base
    if (c.nome == "dois_blocos") num.nx = 20, num.ny = 40;  // 0,25 m x 0,25 m
    num.verbose = 0;
    return num;
}

/// cada combinacao (caso, malha, curva, estado inicial) roda uma vez
const Resultado &Roda(const Caso &c, const Numerica &num) {
    static std::map<std::string, Resultado> guardados;
    std::ostringstream k;
    k << c.nome << ' ' << num.nx << ' ' << num.ny << ' ' << num.rx << ' ' << num.ret << ' ' << num.aG << ' ' << num.bVG
      << ' ' << num.pot << ' ' << num.ci;
    auto it = guardados.find(k.str());
    if (it == guardados.end()) {
        it = guardados.emplace(k.str(), Executa(c, num)).first;
        const Resultado &r = it->second;
        std::cout << "[V3] " << k.str() << ": " << r.passos << " passos, " << r.newton << " Newton, " << r.seg
                  << " s; saida entre y = " << r.ySaida << " e " << r.yAcima << " m; q_esq = " << r.ultimo.qEsq
                  << " m3/s/m (Dupuit-Charny " << r.qDC << ", com a franja " << r.qDCcorr << ")" << std::endl;
    }
    return it->second;
}

std::vector<std::array<REAL, 2>> Ref(const Caso &c) { return LeLinha(std::string(VERIF_DIR) + "/dados/" + c.csvLinha); }

/// espacamento dos nos da face a jusante perto do ponto de saida
REAL Espacamento(const Resultado &r) { return r.yAcima - r.ySaida; }

/// verificacoes comuns a todo regime permanente: convergiu e conserva a massa
void Permanente(const Resultado &r, REAL tolEst) {
    INFO(r.msg);
    REQUIRE(r.ok);
    REQUIRE(r.permanente);
    const REAL qSai = r.ultimo.qJus + r.ultimo.qFace + r.ultimo.qTopo;
    CHECK(std::fabs(r.ultimo.qEsq - qSai) / r.ultimo.qEsq < tolEst);
    // nenhuma agua sai pelo topo no regime permanente (a linha freatica fica abaixo da crista)
    CHECK(std::fabs(r.ultimo.qTopo) < 1.e-6 * r.ultimo.qEsq);
}

} // namespace

// ------------------------------------------------------------------------------------------------- componentes
TEST_CASE("V3: malha da barragem, linha freatica de um campo linear e leitura dos dados", "[V3]") {
    // a malha do quadrilatero: area e comprimentos dos contornos exatos (mapa bilinear com lados retos)
    for (const Caso &c : {CasoRetangular(), CasoInclinada(), CasoDoisBlocos()}) {
        Numerica num = NumTeste(c);
        std::unique_ptr<TPZGeoMesh> g(MalhaBarragem(c, num));
        REAL area = 0.;
        std::map<int, REAL> comp;
        for (int64_t i = 0; i < g->NElements(); ++i) {
            TPZGeoEl *gel = g->Element(i);
            const TPZGeoNode &a = *gel->NodePtr(0), &b = *gel->NodePtr(1);
            if (gel->Dimension() == 2) {
                const TPZGeoNode &cc = *gel->NodePtr(2);
                const REAL det = (b.Coord(0) - a.Coord(0)) * (cc.Coord(1) - a.Coord(1)) -
                                 (cc.Coord(0) - a.Coord(0)) * (b.Coord(1) - a.Coord(1));
                CHECK(det > 0.); // orientacao anti-horaria, nenhum triangulo degenerado
                area += 0.5 * det;
            } else
                comp[gel->MaterialId()] += std::hypot(b.Coord(0) - a.Coord(0), b.Coord(1) - a.Coord(1));
        }
        CHECK(area == Approx(0.5 * (c.xB + c.xD) * c.yE).epsilon(1.e-12));
        CHECK(comp[EBase] == Approx(c.xB).epsilon(1.e-12));
        CHECK(comp[ERight] == Approx(std::hypot(c.xB - c.xD, c.yD)).epsilon(1.e-12));
        CHECK(comp[ETop] == Approx(c.xD).epsilon(1.e-12));
        CHECK(comp[ELeft] == Approx(c.yE).epsilon(1.e-12));
    }
    // linha freatica de p_w = gamma_w (h(x) - y) com h linear: P1 e exato, a isolinha p_w = 0 e y = h(x)
    {
        const Caso c = CasoInclinada();
        Numerica num = NumTeste(c);
        std::unique_ptr<TPZGeoMesh> g(MalhaBarragem(c, num));
        Resultado r;
        auto h = [](REAL x) { return 4.6 - 0.3 * x; };
        for (int64_t i = 0; i < g->NNodes(); ++i) {
            const REAL x = g->NodeVec()[i].Coord(0), y = g->NodeVec()[i].Coord(1);
            r.nosPw.push_back({x, y, gw * (h(x) - y)});
        }
        for (int64_t i = 0; i < g->NElements(); ++i)
            if (g->Element(i)->Dimension() == 2)
                r.tri.push_back({g->Element(i)->NodeIndex(0), g->Element(i)->NodeIndex(1), g->Element(i)->NodeIndex(2)});
        for (REAL x : {0., 0.3, 1., 1.77, 2.5}) CHECK(r.Freatica(x) == Approx(h(x)).margin(1.e-9));
    }
    // os dados digitalizados: as linhas comecam no nivel de montante e terminam no ponto de saida do texto
    const auto l7 = Ref(CasoRetangular()), l14 = Ref(CasoInclinada());
    REQUIRE(l7.size() > 20);
    REQUIRE(l14.size() > 10);
    CHECK(l7.front()[1] == Approx(24.));
    CHECK(l7.back()[0] == Approx(16.));
    CHECK(l7.back()[1] == Approx(12.9));
    CHECK(l14.back()[1] == Approx(3.4));
    CHECK(Interpola(l7, 8.) == Approx(20.50));
}

// ------------------------------------------------------------------------------------------ 4.1 barragem retangular
TEST_CASE("V3: barragem retangular (Sec. 4.1) - vazao de Dupuit-Charny", "[V3]") {
    const Caso c = CasoRetangular();
    const Resultado &r = Roda(c, NumTeste(c));
    Permanente(r, 1.e-4);
    // q = K (H1^2 - H2^2) / (2 L) = 17,5 K = 1,717e-2 m3/s por metro e exato para a superficie livre nitida (prova de
    // Charny). O nosso solo nao saturado conduz um pouco de agua acima da linha freatica: a mesma identidade com o
    // potencial de Kirchhoff (qDCcorr, calculado com as succoes da face a jusante) da a vazao exata do problema de
    // Richards. Tolerancias: contra Dupuit-Charny, 0,3% (discretizacao) + o limite a priori do acrescimo da franja,
    // qFranjaMax / qDC = (topo - H2) |Phi~_inf| / (R qDC) (0,57% com a = 1/kPa: franja de ~0,08 m de agua numa face
    // de 20 m); contra a identidade corrigida, 0,2% (so o erro das reacoes nodais P1 numa malha de 1 m).
    INFO("q_esq = " << r.ultimo.qEsq << " q_DC = " << r.qDC << " q_DC_corr = " << r.qDCcorr
                    << " franja max = " << r.qFranjaMax);
    CHECK(r.qDC == Approx(17.5 * c.K).epsilon(1.e-12));
    CHECK(std::fabs(r.ultimo.qEsq - r.qDC) / r.qDC < 0.003 + r.qFranjaMax / r.qDC);
    CHECK(r.ultimo.qEsq >= r.qDC); // a franja capilar so acrescenta vazao
    CHECK(std::fabs(r.ultimo.qEsq - r.qDCcorr) / r.qDCcorr < 0.002);
}

TEST_CASE("V3: barragem retangular (Sec. 4.1) - ponto de saida e linha freatica", "[V3]") {
    const Caso c = CasoRetangular();
    const Resultado &r = Roda(c, NumTeste(c));
    Permanente(r, 1.e-4);
    // ponto de saida: o artigo da y_s = 8,9 m acima do nivel de jusante (y = 12,9 m). O nosso esta entre o no mais
    // alto da face de percolacao e o no seguinte (sem fluxo, p_w < 0); estimativa = ponto medio. Tolerancia = meio
    // espacamento dos nos da face (resolucao) + 0,3 m (diferenca entre as solucoes publicadas na Fig. 7 junto a face e
    // a franja capilar do nosso modelo)
    const REAL dy = Espacamento(r);
    INFO("saida entre " << r.ySaida << " e " << r.yAcima << " m; artigo " << c.ySaidaArtigo << " m");
    CHECK(dy > 0.);
    CHECK(std::fabs(r.ySaidaEst() - c.ySaidaArtigo) <= 0.5 * dy + 0.3);
    // a face de percolacao comeca no nivel de jusante e e continua ate a saida
    for (const auto &f : r.face) {
        if (f.y > c.H2 + 1.e-9 && f.y <= r.ySaida + 1.e-9) {
            CHECK(f.tipo == 1);
            CHECK(f.q > 0.);
        }
        if (f.y > r.ySaida + 1.e-9) CHECK(f.pw < 0.);
    }
    // linha freatica contra a Fig. 7 digitalizada (curva "Present"): tolerancia de 0,35 m = leitura da figura (0,05 m)
    // + dispersao das solucoes publicadas (~0,2 m) + malha de 1 m com a franja capilar (~0,1 m)
    const auto ref = Ref(c);
    for (REAL x : {2., 4., 6., 8., 10., 12., 14., 15.}) {
        const REAL y = r.Freatica(x), ya = Interpola(ref, x);
        INFO("x = " << x << ": nosso " << y << " artigo " << ya);
        CHECK(std::fabs(y - ya) <= 0.35);
    }
}

TEST_CASE("V3: barragem retangular (Sec. 4.1) - sensibilidade a curva de retencao", "[V3]") {
    // o modelo nao saturado tende a superficie livre nitida quando a curva fica mais ingreme (a maior): a vazao a mais
    // da franja capilar cai com 1/a e o ponto de saida quase nao muda. Tres curvas na malha de 1 m: a = 0,5, 1 e 2/kPa
    const Caso c = CasoRetangular();
    std::vector<REAL> as = {0.5, 1., 2.};
    std::vector<REAL> excesso, ys;
    for (REAL a : as) {
        Numerica num = NumTeste(c);
        num.aG = a;
        const Resultado &r = Roda(c, num);
        Permanente(r, 1.e-4);
        excesso.push_back(r.ultimo.qEsq / r.qDC - 1.);
        ys.push_back(r.ySaidaEst());
        // a identidade de Charny corrigida vale para qualquer curva (erro de discretizacao)
        CHECK(std::fabs(r.ultimo.qEsq - r.qDCcorr) / r.qDCcorr < 0.002);
        std::cout << "[V3] sensibilidade: a = " << a << " 1/kPa: q / q_DC - 1 = " << excesso.back()
                  << ", saida (ponto medio) " << ys.back() << " m" << std::endl;
    }
    // o excesso de vazao diminui quando a curva fica mais ingreme (a franja ~ 1/a) e fica abaixo de 1% para
    // a >= 0,5/kPa (limite a priori da franja com a = 0,5/kPa: 1,1%)
    CHECK(excesso[0] > excesso[1]);
    CHECK(excesso[1] > excesso[2]);
    for (REAL e : excesso) CHECK(e < 0.01);
    // o ponto de saida varia menos que um espacamento dos nos da face (1 m) entre as tres curvas
    CHECK(*std::max_element(ys.begin(), ys.end()) - *std::min_element(ys.begin(), ys.end()) <= 1. + 1.e-9);
}

// ---------------------------------------------------------------------------------- 4.2 face a jusante inclinada
TEST_CASE("V3: barragem de face inclinada (Sec. 4.2) - face de percolacao de 2,4 m", "[V3]") {
    const Caso c = CasoInclinada();
    const Resultado &r = Roda(c, NumTeste(c));
    Permanente(r, 1.e-4);
    // saida em y = 3,4 m (face de percolacao de 2,4 m acima do nivel de jusante de 1 m). Tolerancia = meio
    // espacamento dos nos da face (0,125 m nesta malha) + 0,15 m (na Fig. 14 a linha encontra a face perto de
    // y = 3,3 m: a propria figura difere ~0,1 m do texto)
    const REAL dy = Espacamento(r);
    INFO("saida entre " << r.ySaida << " e " << r.yAcima << " m; artigo " << c.ySaidaArtigo << " m");
    CHECK(std::fabs(r.ySaidaEst() - c.ySaidaArtigo) <= 0.5 * dy + 0.15);
    CHECK(std::fabs((r.ySaidaEst() - c.H2) - 2.4) <= 0.5 * dy + 0.15);
    // linha freatica contra a Fig. 14 digitalizada: 0,12 m = leitura (0,05 m) + franja e malha (~0,07 m)
    const auto ref = Ref(c);
    for (REAL x : {0.5, 1., 1.5, 2., 2.5, 3.}) {
        const REAL y = r.Freatica(x), ya = Interpola(ref, x);
        INFO("x = " << x << ": nosso " << y << " artigo " << ya);
        CHECK(std::fabs(y - ya) <= 0.12);
    }
}

// --------------------------------------------------------------------------------------------- casos demorados
TEST_CASE("V3: independencia do estado inicial (modelos A e linear)", "[.lento][V3]") {
    // como na Fig. 9 do artigo: estados iniciais diferentes levam ao mesmo regime permanente
    const Caso c = CasoRetangular();
    Numerica a = NumTeste(c), b = NumTeste(c);
    b.ci = "linear";
    const Resultado &ra = Roda(c, a), &rb = Roda(c, b);
    Permanente(ra, 1.e-4);
    Permanente(rb, 1.e-4);
    CHECK(rb.ySaida == Approx(ra.ySaida).margin(1.e-9));
    CHECK(rb.ultimo.qEsq == Approx(ra.ultimo.qEsq).epsilon(1.e-3));
    for (REAL x : {4., 8., 12., 15.}) CHECK(std::fabs(ra.Freatica(x) - rb.Freatica(x)) < 0.02);
}

TEST_CASE("V3: malhas do caso completo (retangular 0,5 m e inclinada 0,125 m)", "[.lento][V3]") {
    {
        const Caso c = CasoRetangular();
        Numerica num = NumTeste(c);
        num.nx = 32, num.ny = 48;
        const Resultado &r = Roda(c, num);
        Permanente(r, 1.e-4);
        INFO("q_esq = " << r.ultimo.qEsq << " q_DC = " << r.qDC << " q_DC_corr = " << r.qDCcorr << "; saida entre "
                        << r.ySaida << " e " << r.yAcima);
        CHECK(std::fabs(r.ultimo.qEsq - r.qDC) / r.qDC < 0.003 + r.qFranjaMax / r.qDC);
        CHECK(std::fabs(r.ultimo.qEsq - r.qDCcorr) / r.qDCcorr < 0.001);
        CHECK(std::fabs(r.ySaidaEst() - c.ySaidaArtigo) <= 0.5 * Espacamento(r) + 0.3);
        const auto ref = Ref(c);
        for (REAL x : {2., 4., 6., 8., 10., 12., 14., 15., 15.5}) CHECK(std::fabs(r.Freatica(x) - Interpola(ref, x)) <= 0.3);
    }
    {
        const Caso c = CasoInclinada();
        Numerica num = NumTeste(c);
        num.nx = 28, num.ny = 40;
        const Resultado &r = Roda(c, num);
        Permanente(r, 1.e-4);
        INFO("saida entre " << r.ySaida << " e " << r.yAcima);
        CHECK(std::fabs(r.ySaidaEst() - c.ySaidaArtigo) <= 0.5 * Espacamento(r) + 0.15);
        const auto ref = Ref(c);
        for (REAL x : {0.5, 1., 1.5, 2., 2.5, 3., 3.2}) CHECK(std::fabs(r.Freatica(x) - Interpola(ref, x)) <= 0.1);
    }
}

TEST_CASE("V3: barragem em dois blocos (Sec. 4.3)", "[.lento][V3]") {
    // bloco esquerdo 10 vezes menos permeavel: a linha freatica cai quase na vertical ao longo da interface. Artigo:
    // cruza a interface (x = 2,5 m) a 8,05 m e sai 0,46 m acima do nivel de jusante (y = 2,46 m). Tolerancias: 0,3 m na
    // interface (a queda e quase vertical: um erro pequeno em x vira um erro grande em y; varias referencias da
    // Fig. 17 diferem entre si) e meio espacamento + 0,15 m na saida. A cota na interface e a linha do bloco esquerdo
    // extrapolada das suas duas ultimas colunas de nos (yInterface): na propria interface o nosso p_w fica na franja
    // capilar (entre -1 e 0 kPa) no trecho em que a linha desce por ela, e a isolinha p_w = 0 marca o pe desse trecho
    const Caso c = CasoDoisBlocos();
    const Resultado &r = Roda(c, NumTeste(c));
    Permanente(r, 1.e-4);
    INFO("interface: nosso " << r.yInterface << " (p_w = 0 na interface: " << r.yInterface0 << ") artigo "
                             << c.yInterfaceArtigo << "; saida entre " << r.ySaida
                             << " e " << r.yAcima << " artigo " << c.ySaidaArtigo);
    CHECK(std::fabs(r.yInterface - c.yInterfaceArtigo) <= 0.3);
    CHECK(std::fabs(r.ySaidaEst() - c.ySaidaArtigo) <= 0.5 * Espacamento(r) + 0.15);
    // a identidade de Charny vale para blocos em serie: q = (H1^2 - H2^2) / (2 sum(L_i / K_i)) = 1,712e-3 m3/s por
    // metro; mesmas tolerancias do caso 4.1 (0,3% + o limite a priori da franja, que aqui e ~1,3%: a franja na face
    // de 8 m do bloco permeavel pesa mais contra a vazao pequena que o bloco esquerdo deixa passar; 0,2% contra a
    // formula com a franja)
    INFO("q_esq = " << r.ultimo.qEsq << " q_DC = " << r.qDC << " q_DC_corr = " << r.qDCcorr
                    << " franja max = " << r.qFranjaMax);
    CHECK(r.qDC == Approx(96. / (2. * (2.5 / (0.1 * c.K) + 2.5 / c.K))).epsilon(1.e-12));
    CHECK(std::fabs(r.ultimo.qEsq - r.qDC) / r.qDC < 0.003 + r.qFranjaMax / r.qDC);
    CHECK(std::fabs(r.ultimo.qEsq - r.qDCcorr) / r.qDCcorr < 0.002);
}
