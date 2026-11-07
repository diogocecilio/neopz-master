// V2_TaludesCho: testes Catch2 (os demorados tem a marca [.lento]). Valores publicados: Liao et al. (2025) e Cho (2010);
// a tabela "artigo x nosso resultado" e a discussao estao no README.md.
#define CATCH_CONFIG_MAIN
#include <catch2/catch.hpp>

#include "problema.h"

#include <fstream>
#include <map>
#include <unistd.h>

namespace {

v2::Config CfgTeste(const std::string &caso, double h) {
    v2::Config cf = v2::Config::Padrao(caso);
    cf.h = h;
    return cf;
}

// FS deterministico (as medias) de uma configuracao
double FSDeterministico(const v2::Config &cf, double *seg = nullptr) {
    v2::Simulador sim(cf);
    const auto r = sim.Roda(0, false);
    REQUIRE(r.situacao == "ok");
    if (seg) *seg = r.seg;
    return r.srm.FS;
}

// O caso deterministico do SRM por bissecao com a malha de Liao (h = 0,5 m), calculado uma vez por processo: os testes
// do valor publicado e os do arc-length comparam com ele
const v2::Realizacao &DeterministicoBissecao(const std::string &caso) {
    static std::map<std::string, v2::Realizacao> cache;
    auto it = cache.find(caso);
    if (it == cache.end()) {
        v2::Simulador sim(CfgTeste(caso, 0.5));
        it = cache.emplace(caso, sim.Roda(0, false)).first;
    }
    REQUIRE(it->second.situacao == "ok");
    return it->second;
}

// dados/autovalores_liao.csv: caso -> autovalores e epsilon (digitalizados das Figs. 3 e 8)
std::map<std::string, std::vector<std::pair<double, double>>> AutovaloresLiao() {
    std::map<std::string, std::vector<std::pair<double, double>>> m;
    std::ifstream f(std::string(VERIF_DIR) + "/dados/autovalores_liao.csv");
    REQUIRE(f.good());
    std::string l;
    while (std::getline(f, l)) {
        if (l.empty() || l[0] == '#' || l.rfind("caso", 0) == 0) continue;
        std::stringstream ss(l);
        std::string caso, M, lam, eps;
        std::getline(ss, caso, ',');
        std::getline(ss, M, ',');
        std::getline(ss, lam, ',');
        std::getline(ss, eps, ',');
        m[caso].push_back({std::atof(lam.c_str()), eps.empty() ? -1. : std::atof(eps.c_str())});
    }
    return m;
}

} // namespace

// ---------------------------------------------------------------------------------------------------- KL (sem MEF)

// O Nystrom (Gauss-Legendre + Jacobi) e conferido no nucleo exponencial, que tem autopares em forma fechada
// (Comum/KLSeparavel.h): o nucleo tem bico em s = t, entao a quadratura converge como O(n^-2); com n = 200 o erro
// relativo dos 6 primeiros autovalores fica abaixo de 1e-3 (medido: 8e-5 no primeiro, 9e-4 no sexto) e o das
// autofuncoes abaixo de 1e-3 (a menos do sinal, arbitrario). No nucleo gaussiano (analitico) a convergencia e
// exponencial: 80 e 160 pontos dao os mesmos autovalores a 1e-9 (relativo ao primeiro). O traco discreto e o
// comprimento do intervalo (sum lambda = D, a variancia unitaria integrada).
TEST_CASE("V2: Nystrom reproduz os autopares do nucleo exponencial e converge no gaussiano", "[V2]") {
    const double D = 10., L = 2.5;
    const auto bn = v2kl::BaseNystrom(D, [L](double r) { return std::exp(-std::fabs(r) / L); }, 200);
    const auto be = v2kl::BaseExponencial(D, L, 10);
    for (int i = 0; i < 6; ++i) {
        CHECK(bn.lambda[i] == Approx(be.lambda[i]).epsilon(1.e-3));
        for (double s : {-4.9, -2.2, 0.3, 1.3, 4.4}) CHECK(std::fabs(bn.f(i, s)) == Approx(std::fabs(be.f(i, s))).margin(1.e-3));
    }
    double traco = 0.;
    for (double l : bn.lambda) traco += l;
    CHECK(traco == Approx(D).epsilon(1.e-6));
    const auto g1 = v2kl::BaseGaussiana(D, L, 80), g2 = v2kl::BaseGaussiana(D, L, 160);
    for (int i = 0; i < 12; ++i) CHECK(std::fabs(g1.lambda[i] - g2.lambda[i]) < 1.e-9 * g1.lambda[0]);
}

// Liao et al. (2025), Figs. 3 e 8 (dados/autovalores_liao.csv, digitalizados): autovalores da KL de c_u (cho1:
// h_x = 25 m, h_z = 2,5 m) e de c, phi (cho2: 20 m e 2 m) e o numero de termos para epsilon >= 95% (M = 9 e 18, secs.
// 4.1 e 4.2). Com o nucleo gaussiano (rho2) no retangulo que contem o talude (30 x 10 e 30 x 15) obtemos M = 9 e 18 e
// lambda_1 a menos de 3% (os outros ficam ate ~17% abaixo: a fracao epsilon de Liao implica um traco de ~317 e ~468,
// nao 300 e 450). Com o retangulo aumentado em meia celula da grade de 0,5 m (KL_MARGEM = 0,25: 30,5 x 10,5 e
// 30,5 x 15,5, a discretizacao nodal de 61 x 21 e 61 x 31 nos com 0,25 m2 cada) os autovalores coincidem com os
// digitalizados a menos da leitura (+-0,5: ~1% nos primeiros, ate ~10% nos menores), tolerancia de 4% + 0,5, e
// epsilon a 0,02. O nucleo exponencial (rho1) com as mesmas distancias da lambda_1 = 81 e 70 (contra 100 e 78) e
// M = 133 e 295: as Figs. 3 e 8 sao do nucleo gaussiano.
TEST_CASE("V2: autovalores da KL gaussiana contra as Figs. 3 e 8 de Liao", "[V2]") {
    const auto dig = AutovaloresLiao();
    REQUIRE(dig.count("cho1"));
    REQUIRE(dig.count("cho2"));
    for (const std::string caso : {"cho1", "cho2"}) {
        const v2::Caso cs = v2::CasoPorNome(caso);
        const int Mliao = caso == "cho1" ? 9 : 18;
        const auto &d = dig.at(caso);
        const auto kl = v2kl::CriaKL("gauss", 0., 0., cs.L, cs.zc, cs.hx, cs.hz, 0.95, 2000);
        CHECK(kl.M() == Mliao);
        CHECK(kl.Lambda(0) == Approx(d[0].first).epsilon(0.03));
        // o retangulo da discretizacao nodal de Liao (50 termos, para comparar todos os digitalizados)
        const double m = 0.25, Dx = cs.L + 2 * m, Dy = cs.zc + 2 * m;
        CHECK(v2kl::CriaKL("gauss", -m, -m, Dx, Dy, cs.hx, cs.hz, 0.95, 2000).M() == Mliao);
        const auto kn = v2kl::CriaKL("gauss", -m, -m, Dx, Dy, cs.hx, cs.hz, 1., 50);
        REQUIRE(kn.M() >= (int)d.size());
        double soma = 0.;
        for (size_t k = 0; k < d.size(); ++k) {
            CHECK(kn.Lambda((int)k) == Approx(d[k].first).epsilon(0.04).margin(0.5));
            soma += kn.Lambda((int)k);
            if (d[k].second > 0.) CHECK(soma / (Dx * Dy) == Approx(d[k].second).margin(0.02));
        }
        // nucleo exponencial: outra curva
        const auto ke = v2kl::CriaKL("exp", 0., 0., cs.L, cs.zc, cs.hx, cs.hz, 0.95, 2000);
        CHECK(ke.M() > 5 * Mliao);
        CHECK(ke.Lambda(0) < 0.95 * d[0].first);
    }
}

// Estatistica dos campos gerados (sem MEF), para os dois nucleos, em pontos com defasagens horizontais (5 e 10 m) e
// verticais (0,5, 1 e 2,5 m) a partir de (15, 5), com N = 20 000 realizacoes. Tolerancias: (a) contra a propria
// expansao truncada (media de ln c = mu_ln, desvio de ln c = sigma_ln sqrt(var_M), covariancia C_M): 5 erros-padrao
// do estimador (media: sigma/sqrt(N); desvio: sigma/sqrt(2N); covariancia: (1 + rho^2)^0,5 / sqrt(N)); (b) contra o
// nucleo exato rho(dx, dz): o resto da expansao e um campo de variancia 1 - var_M em cada ponto, logo, por
// Cauchy-Schwarz, |rho - C_M| <= sqrt((1 - var_M(p))(1 - var_M(q))), mais os 5 erros-padrao; (c) CoV de c contra o
// alvo 0,3: entre o valor truncado sqrt(exp(sigma_ln^2 var_M) - 1) e 0,3, mais 5 erros-padrao. No cho2, a correlacao
// cruzada de ln c e ln phi no mesmo ponto e r = -0,5 exatamente (mesma base, eq. 12): 5 erros-padrao,
// (1 - r^2)/sqrt(N).
TEST_CASE("V2: campos lognormais por KL: media, CoV e correlacao (nucleos gauss e exp)", "[V2]") {
    const std::vector<double> px = {15., 20., 25., 15., 15., 15.}, py = {5., 5., 5., 5.5, 6., 7.5};
    const int N = 20000;
    for (const std::string caso : {"cho1", "cho2"})
        for (const std::string nuc : {"gauss", "exp"}) {
            INFO(caso << " " << nuc);
            const v2::Caso cs = v2::CasoPorNome(caso);
            const v2::GeradorCampos g(cs, nuc, 0.95, 2000, 12345ULL, px, py);
            const int np = (int)px.size();
            std::vector<double> s1(np, 0.), s2(np, 0.), sc(np, 0.), sc2(np, 0.), sx(np, 0.);
            double sxy = 0., sy2 = 0.;
            std::vector<double> zc, zp;
            for (int r = 0; r < N; ++r) {
                g.Gaussianos(r, zc, zp);
                for (int p = 0; p < np; ++p) {
                    s1[p] += zc[p], s2[p] += zc[p] * zc[p], sx[p] += zc[0] * zc[p];
                    const double c = klsep::Lognormal(cs.c, cs.covC)(zc[p]);
                    sc[p] += c, sc2[p] += c * c;
                }
                if (cs.phiAleatorio) sxy += zc[0] * zp[0], sy2 += zp[0] * zp[0];
            }
            const klsep::Lognormal ln(cs.c, cs.covC);
            for (int p = 0; p < np; ++p) {
                const double vM = g.P.var[p], m = s1[p] / N, v = s2[p] / N - m * m;
                // (a) ln c = mu_ln + sigma_ln Z: media e desvio de Z
                CHECK(std::fabs(m) < 5. * std::sqrt(vM / N));
                CHECK(std::fabs(std::sqrt(v) - std::sqrt(vM)) < 5. * std::sqrt(vM / (2. * N)));
                // (c) CoV de c
                const double mc = sc[p] / N, cov = std::sqrt(std::max(0., sc2[p] / N - mc * mc)) / mc;
                const double covM = std::sqrt(std::exp(ln.sigma_ln * ln.sigma_ln * vM) - 1.);
                const double se = cs.covC * 1.5 / std::sqrt(N); // erro-padrao aproximado do CoV de uma lognormal
                CHECK(cov > covM - 5. * se);
                CHECK(cov < cs.covC + 5. * se);
                CHECK(mc == Approx(std::exp(ln.mu_ln + 0.5 * ln.sigma_ln * ln.sigma_ln * vM)).epsilon(5. * cs.covC / std::sqrt(N)));
                // (a) e (b) covariancia com o ponto 0
                const double C = sx[p] / N - (s1[0] / N) * m, CM = g.kl->Covariancia(g.P, 0, p);
                const double rho = v2kl::NucleoExato(nuc, px[p] - px[0], py[p] - py[0], cs.hx, cs.hz);
                const double se2 = std::sqrt((1. + rho * rho) / N);
                CHECK(std::fabs(C - CM) < 5. * se2);
                CHECK(std::fabs(C - rho) < std::sqrt((1. - g.P.var[0]) * (1. - vM)) + 5. * se2);
            }
            std::ostringstream o;
            o << caso << " " << nuc << ": M = " << g.M() << ", var_M(15,5) = " << g.P.var[0] << "; c: media "
              << sc[0] / N << ", CoV " << std::sqrt(sc2[0] / N - (sc[0] / N) * (sc[0] / N)) / (sc[0] / N)
              << "; correlacao de ln c (amostra / KL truncada / nucleo) com dx = 5, 10 e dz = 0,5, 1, 2,5:";
            for (int p = 1; p < np; ++p) {
                const double m = s1[p] / N, C = sx[p] / N - (s1[0] / N) * m;
                const double v0 = s2[0] / N - (s1[0] / N) * (s1[0] / N), vp = s2[p] / N - m * m;
                o << " " << std::setprecision(3) << C / std::sqrt(v0 * vp) << "/"
                  << g.kl->Covariancia(g.P, 0, p) / std::sqrt(g.P.var[0] * g.P.var[p]) << "/"
                  << v2kl::NucleoExato(nuc, px[p] - px[0], py[p] - py[0], cs.hx, cs.hz);
            }
            if (cs.phiAleatorio) {
                const double r = (sxy / N) / std::sqrt((s2[0] / N) * (sy2 / N));
                CHECK(r == Approx(cs.rho).margin(5. * (1. - cs.rho * cs.rho) / std::sqrt(N)));
                o << "; correlacao cruzada ln c x ln phi = " << r;
            }
            WARN(o.str());
        }
}

// ------------------------------------------------------------------------------------------------------------ malha

// A geometria reconstruida das Figs. 2 e 6: com h = 0,5 m a malha tem os nos da grade de Liao (981 no cho1 e 1281 no
// cho2, os numeros das secoes 4.1 e 4.2; exato: e uma contagem), os triangulos cobrem a secao (225 e 300 m2, a 1e-9)
// e a soma dos pesos de integracao dos pontos (o que pondera a estatistica dos campos) da a mesma area.
TEST_CASE("V2: malha: os nos da grade de Liao e a area da secao", "[V2]") {
    for (const std::string caso : {"cho1", "cho2"}) {
        const v2::Caso cs = v2::CasoPorNome(caso);
        v2::Modelo m(cs, 0.5);
        CHECK(m.gmesh->NNodes() == (caso == "cho1" ? 981 : 1281));
        double A = 0.;
        for (TPZGeoEl *gel : m.gmesh->ElementVec())
            if (gel && gel->Dimension() == 2) A += gel->Volume();
        CHECK(A == Approx(cs.AreaSecao()).epsilon(1.e-9));
        CHECK(m.AreaPontos() == Approx(cs.AreaSecao()).epsilon(1.e-9));
        CHECK(m.ips.size() == 3 * (size_t)m.NElementos());
    }
}

// ---------------------------------------------------------------------------------------- FS deterministico (MEF)

// cho1 (Liao, sec. 4.1): FS = 1,369 (SIGMA/W + SLOPE/W) e 1,366 (Jiang et al. 2014), com c_u = 23 kPa e a malha de
// 0,5 m de Liao. O nosso FS vem do SRM (c_u/F ate o Newton nao convergir), outro metodo: o equilibrio limite nos
// circulos da mesma geometria (scripts/equilibrio_limite.py, exato para phi = 0) da 1,354 (circulo tangente a base
// rigida), e o SRM pode achar um mecanismo nao circular um pouco mais critico. Tolerancia de 3% em torno de 1,369:
// ~1% da diferenca de metodo (LEM circular x SRM), ~1% da discretizacao (de h = 1 para 0,5 m o FS caiu 0,9%; ver
// README) e ~0,7% da bissecao (FS_TOL = 0,01). Contra o nosso LEM: 2%.
TEST_CASE("V2: cho1 deterministico: FS por SRM contra 1,369 (Liao) e 1,366 (Jiang)", "[V2]") {
    const auto &r = DeterministicoBissecao("cho1");
    const double FS = r.srm.FS, seg = r.seg;
    WARN("cho1, h = 0,5 m: FS = " << FS << " (" << seg << " s); Liao 1,369, Jiang 1,366, LEM circular 1,354");
    CHECK(FS == Approx(1.369).epsilon(0.03));
    CHECK(FS == Approx(1.366).epsilon(0.03));
    CHECK(FS == Approx(1.354).epsilon(0.02));
}

// cho2 (Liao, sec. 4.2): FS = 1,205 (Bishop; 1,204 em Cho 2010) e 1,284 (elementos finitos, metodo nao informado).
// O SRM com fluxo associado, num talude homogeneo, fica perto do equilibrio limite (Griffiths e Lane 1999): o nosso
// Bishop nesta geometria da 1,215 (scripts/equilibrio_limite.py). Tolerancia de 3% em torno de 1,205 (as mesmas
// parcelas do cho1). O valor de 1,284 nao e alcancado: o SRM fica ~5% abaixo dele (ver a discussao no README); o
// teste so confere que o nosso FS esta entre o Bishop menos 3% e esse valor.
TEST_CASE("V2: cho2 deterministico: FS por SRM contra 1,205 (Bishop) e 1,284 (MEF)", "[V2]") {
    const auto &r = DeterministicoBissecao("cho2");
    const double FS = r.srm.FS, seg = r.seg;
    WARN("cho2, h = 0,5 m: FS = " << FS << " (" << seg << " s); Liao Bishop 1,205, MEF 1,284; nosso Bishop 1,215");
    CHECK(FS == Approx(1.205).epsilon(0.03));
    CHECK(FS == Approx(1.215).epsilon(0.02));
    CHECK(FS > 1.205 * 0.97);
    CHECK(FS < 1.284);
}

// ------------------------------------------------------------------------------- FS por arc-length (SRM e GIM)

// O SRM e o GIM por arc-length (FactorSRMArc e FactorGIM de Comum/RahardjoDriver.h, o codigo do LearningMohrCoulomb2
// do branch UnsaturatedSoil, commit 2f92365, com os seus parametros padrao) no cho1, malha de Liao (h = 0,5 m):
// (1) SRM contra 1,369 e 1,366: 3%, as parcelas do teste da bissecao (metodo e malha; a da parada agora e 0,1%).
// (2) SRM contra a bissecao na mesma malha, que tem o mesmo limite discreto F*: a bissecao da o ultimo F convergido,
//     em [F* - 0,01, F*] se as falhas sao genuinas (FS_TOL = 0,01); o arco da o maior F do caminho, abaixo de F* pelo
//     erro da parada (dF < 0,001 por incremento; no estudo de 43 casos de 2f92365, -0,46% no pior caso). Logo
//     FS_arco - FS_bissecao fica entre -0,5% FS e 0,01 + 0,5% FS.
// (3) Sem atrito, GIM e SRM tem o mesmo limite, tambem no problema discreto: se sigma equilibra o peso lambda gamma com
//     c_u, sigma / lambda equilibra gamma com c_u / lambda (o criterio de Tresca e homogeneo de grau 1 em sigma e c_u, e
//     o limite de um material perfeitamente plastico associado nao depende do caminho). A diferenca vem das duas paradas:
//     1% (0,46% + 0,41%, os piores erros do SRM e do GIM no estudo de 2f92365).
TEST_CASE("V2: cho1: SRM e GIM por arc-length contra 1,369, a bissecao e entre si (phi = 0)", "[V2]") {
    const double FSb = DeterministicoBissecao("cho1").srm.FS;
    v2::Simulador sim(CfgTeste("cho1", 0.5));
    sim.mod->PropriedadesMedias();
    const auto srm = v2::FatorArco(*sim.mod->an, false);
    const auto gim = v2::FatorArco(*sim.mod->an, true);
    WARN("cho1, h = 0,5 m: SRM por arco " << srm.FS << " (" << srm.saida.incrementos << " incrementos, " << srm.saida.cortes
         << " cortes, " << srm.saida.iteracoes << " iteracoes, " << srm.saida.motivo << ", " << srm.seg << " s); GIM por arco "
         << gim.FS << " (" << gim.saida.incrementos << " incrementos, " << gim.saida.cortes << " cortes, "
         << gim.saida.iteracoes << " iteracoes, " << gim.saida.motivo << ", " << gim.seg << " s); bissecao " << FSb << " ("
         << DeterministicoBissecao("cho1").srm.newton << " iteracoes); Liao 1,369");
    CHECK(srm.saida.motivo == "parada");
    CHECK(gim.saida.motivo == "parada");
    CHECK(srm.FS == Approx(1.369).epsilon(0.03));
    CHECK(srm.FS == Approx(1.366).epsilon(0.03));
    CHECK(srm.FS - FSb > -0.005 * FSb);
    CHECK(srm.FS - FSb < 0.01 + 0.005 * FSb);
    CHECK(gim.FS == Approx(srm.FS).epsilon(0.01));
}

// cho2 (c-phi), malha de Liao: (1) SRM por arc-length contra 1,205 (Bishop), 3% como no teste da bissecao, e abaixo de
// 1,284 (MEF); (2) contra a bissecao na mesma malha, com as faixas do cho1. (3) GIM > SRM: o colapso de um talude
// homogeneo so depende de c/(gamma H) e de phi (Taylor): c / (lambda gamma H) = N(phi) no GIM e (c / F) / (gamma H) =
// N(phi_F) no SRM, com tan(phi_F) = tan(phi) / F e N decrescente em phi; com F > 1, phi_F < phi e F < lambda (com
// atrito, aumentar o peso aumenta tambem a resistencia por atrito, e so a coesao resiste ao aumento).
TEST_CASE("V2: cho2: SRM e GIM por arc-length contra 1,205, a bissecao e GIM > SRM (c-phi)", "[V2]") {
    const double FSb = DeterministicoBissecao("cho2").srm.FS;
    v2::Simulador sim(CfgTeste("cho2", 0.5));
    sim.mod->PropriedadesMedias();
    const auto srm = v2::FatorArco(*sim.mod->an, false);
    const auto gim = v2::FatorArco(*sim.mod->an, true);
    WARN("cho2, h = 0,5 m: SRM por arco " << srm.FS << " (" << srm.saida.incrementos << " incrementos, " << srm.saida.cortes
         << " cortes, " << srm.saida.iteracoes << " iteracoes, " << srm.saida.motivo << ", " << srm.seg << " s); GIM por arco "
         << gim.FS << " (" << gim.saida.incrementos << " incrementos, " << gim.saida.cortes << " cortes, "
         << gim.saida.iteracoes << " iteracoes, " << gim.saida.motivo << ", " << gim.seg << " s); bissecao " << FSb << " ("
         << DeterministicoBissecao("cho2").srm.newton << " iteracoes); Liao Bishop 1,205, MEF 1,284");
    CHECK(srm.saida.motivo == "parada");
    CHECK(gim.saida.motivo == "parada");
    CHECK(srm.FS == Approx(1.205).epsilon(0.03));
    CHECK(srm.FS < 1.284);
    CHECK(srm.FS - FSb > -0.005 * FSb);
    CHECK(srm.FS - FSb < 0.01 + 0.005 * FSb);
    CHECK(gim.FS > srm.FS);
}

// phi = 0 (Tresca, o limite do Mohr-Coulomb associado PV2): (1) o FS com phi = 0 e o com phi = 0,001 grau coincidem
// (dentro da tolerancia da bissecao, 0,01): phi = 0 nao e um caso singular do retorno (o apice, que usaria
// cot(phi), nao e alcancado no Tresca); (2) sem atrito o colapso so depende de c_u/(gamma H): FS(2 c_u) = 2 FS(c_u)
// a menos da bissecao (2 x 0,01 / 1,36 + 0,01 / 1,36: 2,2%; usamos 2,5%). Malha de 1 m (o teste e da lei, nao do
// valor do FS).
TEST_CASE("V2: phi = 0 (Tresca) no Mohr-Coulomb associado PV2", "[V2]") {
    auto cf = CfgTeste("cho1", 1.0);
    const double F0 = FSDeterministico(cf);
    cf.phiGraus = 0.001;
    const double Fp = FSDeterministico(cf);
    cf.phiGraus = 0.;
    cf.c = 46.;
    const double F2 = FSDeterministico(cf);
    WARN("cho1, h = 1 m: FS(phi = 0) = " << F0 << ", FS(phi = 0,001 grau) = " << Fp << ", FS(2 c_u) = " << F2);
    REQUIRE(std::isfinite(F0));
    CHECK(F0 > 1.);
    CHECK(std::fabs(Fp - F0) <= 0.01 + 1.e-9);
    CHECK(F2 / F0 == Approx(2.).epsilon(0.025));
}

// ---------------------------------------------------------------------------------------------------- Monte Carlo

// Fumaca do Monte Carlo (malha de 1 m): (1) uma execucao do contrato numa pasta temporaria (CASO=cho1, REAL=1):
// resultado.txt com real=, FS_min=, situacao=ok e seg=; resumo.csv com cabecalho e uma linha (real, FS_min); fs.csv com
// t_h = 0 e FS_SRM; (2) mais duas realizacoes do cho1 no mesmo processo: campos diferentes, FS entre 0,5 e 3 (com
// CoV 0,3 o FS de Liao vai de ~0,7 a ~2,2, Fig. 5b), campos reprodutiveis (a mesma REAL refaz o campo); (3) uma
// realizacao do cho2 com c e phi correlacionados (r = -0,5): phi varia e o SRM termina.
TEST_CASE("V2: Monte Carlo: fumaca (contrato do executavel e realizacoes)", "[V2]") {
    namespace fs = std::filesystem;
    const fs::path pasta = fs::temp_directory_path() / ("v2_fumaca_" + std::to_string(::getpid()));
    fs::create_directories(pasta);
    auto cf = CfgTeste("cho1", 1.0);
    cf.real = 1;
    cf.campos = 1;
    REQUIRE(v2::RodaExecucao(cf, pasta) == 0);
    auto le = [&](const char *a) {
        std::ifstream f(pasta / a);
        std::stringstream ss;
        ss << f.rdbuf();
        return ss.str();
    };
    const std::string res = le("resultado.txt"), resumo = le("resumo.csv"), fsc = le("fs.csv");
    INFO(res);
    CHECK(res.rfind("real=1 ", 0) == 0);
    CHECK(res.find(" FS_min=") != std::string::npos);
    CHECK(res.find(" situacao=ok ") != std::string::npos);
    CHECK(res.find(" seg=") != std::string::npos);
    CHECK(resumo.rfind("real,", 0) == 0);
    CHECK(resumo.find(",FS_min,") != std::string::npos);
    CHECK(std::count(resumo.begin(), resumo.end(), '\n') == 2);
    CHECK(fsc.rfind("t_h,FS_SRM,", 0) == 0);
    CHECK(fsc.find("\n0,") != std::string::npos);
    CHECK(!fs::exists(pasta / "resultado.txt.tmp"));
    const double FS1 = std::atof(res.substr(res.find("FS_min=") + 7).c_str());
    fs::remove_all(pasta);

    v2::Simulador sim(cf);
    std::vector<double> c2a, phi2a, c2b, phi2b;
    sim.Gerador().Gera(2, c2a, phi2a);
    sim.Gerador().Gera(2, c2b, phi2b);
    CHECK(c2a == c2b);
    const auto r2 = sim.Roda(2, true), r3 = sim.Roda(3, true);
    WARN("cho1, h = 1 m, gauss: FS(real 1, 2, 3) = " << FS1 << ", " << r2.srm.FS << ", " << r3.srm.FS << " (" << r2.seg
         << " e " << r3.seg << " s); c medio " << r2.c.media << " e " << r3.c.media << " kPa");
    for (double F : {FS1, r2.srm.FS, r3.srm.FS}) {
        CHECK(F > 0.5);
        CHECK(F < 3.);
    }
    CHECK(r2.c.media != r3.c.media);

    auto c2 = CfgTeste("cho2", 1.0);
    v2::Simulador s2(c2);
    const auto q = s2.Roda(1, true);
    WARN("cho2, h = 1 m, gauss: FS(real 1) = " << q.srm.FS << "; phi de " << q.phi.min << " a " << q.phi.max << " graus");
    CHECK(q.situacao == "ok");
    CHECK(q.phi.max - q.phi.min > 5.);
    CHECK(q.srm.FS > 0.5);
}

// O contrato com FS_METODO=arco (malha de 1 m, cho1, REAL = 1): resultado.txt com metodo=arco e situacao=ok, e o FS da
// mesma realizacao (o mesmo campo) contra o da bissecao com as faixas do teste deterministico do cho1 (-0,5% a
// 0,01 + 0,5%). Com FS_MAX = 1 o caminho termina ao passar de 1: a indicadora FS < 1 e a mesma do caminho completo.
TEST_CASE("V2: Monte Carlo: FS_METODO=arco no contrato e na mesma realizacao da bissecao", "[V2]") {
    namespace fs = std::filesystem;
    const fs::path pasta = fs::temp_directory_path() / ("v2_arco_" + std::to_string(::getpid()));
    fs::create_directories(pasta);
    auto cf = CfgTeste("cho1", 1.0);
    cf.real = 1;
    cf.campos = 1;
    cf.metodo = "arco";
    REQUIRE(v2::RodaExecucao(cf, pasta) == 0);
    std::ifstream f(pasta / "resultado.txt");
    std::stringstream ss;
    ss << f.rdbuf();
    const std::string res = ss.str();
    INFO(res);
    CHECK(res.find(" situacao=ok ") != std::string::npos);
    CHECK(res.find(" metodo=arco") != std::string::npos);
    const double FSa = std::atof(res.substr(res.find("FS_min=") + 7).c_str());
    fs::remove_all(pasta);

    cf.metodo = "bissecao";
    v2::Simulador sim(cf);
    const auto b = sim.Roda(1, true);
    cf.metodo = "arco";
    cf.fsMax = 1.;
    v2::Simulador sa(cf);
    const auto c = sa.Roda(1, true);
    WARN("cho1, h = 1 m, real 1: FS por arco " << FSa << ", por bissecao " << b.srm.FS << " (" << b.srm.newton
         << " iteracoes contra " << c.srm.newton << " do arco censurado em 1, que da " << c.srm.FS << ")");
    CHECK(FSa - b.srm.FS > -0.005 * b.srm.FS);
    CHECK(FSa - b.srm.FS < 0.01 + 0.005 * b.srm.FS);
    CHECK((c.srm.FS < 1.) == (FSa < 1.));
}

// P_f do cho1 contra 10,12% (Liao, Tabela 2: Monte Carlo simples com 10 000 realizacoes do MEF): V2_NMC realizacoes
// (padrao 200) por nucleo, malha V2_HMC (padrao 1 m, para caber num PC; o executavel do cluster usa 0,5 m) e a busca do
// FS censurada em 1 (FS_MAX = 1: a indicadora FS < 1 e a mesma da busca completa). O nucleo gaussiano e o do artigo:
// o teste exige que 10,12% esteja no intervalo de Wilson de 95% da nossa estimativa (com N = 200 e P_f ~ 0,1, cerca de
// +-4 pontos percentuais). O exponencial (rho1) so e informado: com a mesma variancia e campos mais rugosos, a media
// espacial ao longo da superficie de ruptura reduz a dispersao do FS, e a P_f deve ficar abaixo da do gaussiano.
// Grava mc_lento_<nucleo>.csv (real, FS) na pasta corrente.
TEST_CASE("V2: Monte Carlo: P_f do cho1 contra 10,12% (nucleos gauss e exp)", "[V2][.lento]") {
    const int N = v2::EnvI("V2_NMC", 200);
    const double h = v2::EnvD("V2_HMC", 1.0);
    std::map<std::string, std::pair<int, int>> pf;
    for (const std::string nuc : {"gauss", "exp"}) {
        auto cf = CfgTeste("cho1", h);
        cf.nucleo = nuc;
        cf.fsMax = 1.;
        v2::Simulador sim(cf);
        std::ofstream o("mc_lento_" + nuc + ".csv");
        o << "real,FS,c_med,seg\n";
        int falhas = 0;
        const auto t0 = v2::Relogio::now();
        for (int r = 1; r <= N; ++r) {
            const auto q = sim.Roda(r, true);
            REQUIRE(q.situacao == "ok");
            if (q.srm.FS < 1.) ++falhas;
            o << r << ',' << q.srm.FS << ',' << q.c.media << ',' << q.seg << '\n';
            o.flush();
        }
        pf[nuc] = {falhas, N};
        const auto w = v2::Wilson(falhas, N);
        WARN("cho1 " << nuc << " (h = " << h << " m, KL_M = " << sim.Gerador().M() << "): P_f = " << double(falhas) / N
             << " (" << falhas << "/" << N << "), Wilson 95% [" << w.first << ", " << w.second << "], "
             << v2::Segundos(t0) << " s");
    }
    const auto wg = v2::Wilson(pf["gauss"].first, pf["gauss"].second);
    CHECK(0.1012 >= wg.first);
    CHECK(0.1012 <= wg.second);
}
