// V4_TaludeLiu: testes Catch2 contra Liu et al. (2024). Os padrao usam a malha de 1 m (malhas/talude_liu_h10.msh) e
// rodam em cerca de 1 minuto; os demorados (malha de 0,5 m do artigo e um Monte Carlo curto) tem a marca [.lento]:
//   V4_testes                 (padrao)      V4_testes "[.lento]"      (demorados)
#define CATCH_CONFIG_MAIN
#include <catch2/catch.hpp>

#include "problema.h"

#include <cstdlib>

namespace {

/// le um CSV de referencia (linhas com # sao comentarios; a primeira linha restante e o cabecalho) com 2 colunas
std::vector<std::pair<double, double>> LeRef(const std::string &nome) {
    std::ifstream f(std::string(VERIF_DIR) + "/dados/" + nome);
    REQUIRE(f.good());
    std::vector<std::pair<double, double>> v;
    std::string l;
    bool cab = true;
    while (std::getline(f, l)) {
        if (l.empty() || l[0] == '#') continue;
        if (cab) { cab = false; continue; }
        const auto k = l.find(',');
        v.push_back({std::atof(l.substr(0, k).c_str()), std::atof(l.substr(k + 1).c_str())});
    }
    return v;
}

/// cota em que p_w = 0 numa lista (z, p_w) crescente em z (interpolacao linear)
double CotaLencol(const std::vector<std::pair<double, double>> &zp) {
    for (size_t i = 0; i + 1 < zp.size(); ++i)
        if (zp[i].second > 0. && zp[i + 1].second <= 0.)
            return zp[i].first + zp[i].second / (zp[i].second - zp[i + 1].second) * (zp[i + 1].first - zp[i].first);
    return -1.;
}

/// configuracao dos testes: o caso do artigo com a malha pedida, sem mensagens
v4::Config ConfigTeste(const std::string &malha) {
    v4::Config cf;
    cf.malha = malha;
    v4::Verboso() = 0;
    return cf;
}

/// regime permanente e FS deterministicos numa malha, calculados uma vez por execucao dos testes
struct Deterministico {
    v4::InfoPermanente perm;
    v4::InfoFS fs;
    std::vector<std::pair<double, double>> xx, yy; ///< (z, p_w) nas secoes X-X e Y-Y
    std::string pw;                                ///< pw.csv do permanente (para o teste do Monte Carlo)
};
const Deterministico &CasoDeterministico(const std::string &malha) {
    static std::map<std::string, Deterministico> cache;
    auto it = cache.find(malha);
    if (it != cache.end()) return it->second;
    Deterministico d;
    v4::Talude t(ConfigTeste(malha));
    d.perm = t.Permanente();
    for (double z : v4::CotasSecao(v4::kXSecaoX)) d.xx.push_back({z, t.PwEm(v4::kXSecaoX, z)});
    for (double z : v4::CotasSecao(v4::kXSecaoY)) d.yy.push_back({z, t.PwEm(v4::kXSecaoY, z)});
    d.pw = (std::filesystem::temp_directory_path() / ("v4_pw_" + malha + ".csv")).string();
    t.GravaPw(d.pw);
    d.fs = t.FS();
    return cache[malha] = d;
}

/// compara um perfil com a referencia digitalizada; tolSat: abaixo do lencol do artigo, tolNSat: acima
void ComparaPerfil(const std::vector<std::pair<double, double>> &nosso, const std::vector<std::pair<double, double>> &ref,
                   double tolSat, double tolNSat, const char *secao) {
    REQUIRE(nosso.size() >= ref.size());
    double maxSat = 0., maxNSat = 0.;
    for (const auto &r : ref) {
        double pw = std::numeric_limits<double>::quiet_NaN();
        for (const auto &n : nosso) if (std::fabs(n.first - r.first) < 1.e-9) pw = n.second;
        INFO(secao << ": z = " << r.first << " m, artigo " << r.second << " kPa, nosso " << pw << " kPa");
        REQUIRE(std::isfinite(pw));
        const double d = std::fabs(pw - r.second);
        (r.second > 0. ? maxSat : maxNSat) = std::max(r.second > 0. ? maxSat : maxNSat, d);
        CHECK(d <= (r.second > 0. ? tolSat : tolNSat));
    }
    WARN(secao << ": maior diferenca " << maxSat << " kPa abaixo do lencol e " << maxNSat << " kPa acima; lencol em z = "
               << CotaLencol(nosso) << " m (artigo " << CotaLencol(ref) << " m)");
    CHECK(std::fabs(CotaLencol(nosso) - CotaLencol(ref)) <= 0.15);
}

/// FXH (eq. 9 de Liu et al. 2024) com a retencao das eqs. 7-8, integrada pela regra do trapezio em y = ln(psi)
double KrFXH(double psi, double aev) {
    const double a = 5., n = 2., m = 1., psir = 3000.;
    auto th = [&](double s) {
        return (1. - std::log(1. + s / psir) / std::log(1. + 1.e6 / psir)) / std::pow(std::log(M_E + std::pow(s / a, n)), m);
    };
    auto dth = [&](double s) { return (th(s * (1. + 1.e-6)) - th(s * (1. - 1.e-6))) / (2.e-6 * s); };
    auto integral = [&](double lo, double tref) {
        const int N = 20000;
        const double b = std::log(1.e6), h = (b - lo) / N;
        double s = 0.;
        for (int i = 0; i <= N; ++i) {
            const double y = lo + i * h, ey = std::exp(y), f = (th(ey) - tref) / ey * dth(ey);
            s += (i == 0 || i == N ? 0.5 : 1.) * f;
        }
        return s * h;
    };
    if (psi <= aev) return 1.;
    return integral(std::log(psi), th(psi)) / integral(std::log(aev), 1.);
}

} // namespace

// ---------------------------------------------------------------------------------------------------------------------
// Verificacao 1: as funcoes hidraulicas. A retencao da biblioteca (Fredlund-Xing sem C(psi)) contra a eq. 7 do artigo
// (com C(psi), psi_r = 3000 kPa): C >= 0,9989 ate 20 kPa, entao |dS| <= 0,002. A FXH (eq. 9) nao existe na biblioteca:
// k_r = S^p com p = 5,57, o ajuste por minimos quadrados de log10 k_r entre psi_aev = 2,58 kPa (tangente no ponto de
// inflexao) e 20 kPa (scripts/ajuste_kr.py), cujo maior residuo e 9 % (em 20 kPa): tolerancia de 0,05 em log10 k_r
// (12 %). O patamar de infiltracao permanente (gradiente unitario, k_r(psi) = q/k_s = 0,01) do ajuste fica a 2 % do da
// FXH (13,3 x 13,4 kPa).
TEST_CASE("V4: retencao de Fredlund-Xing e k_r = S^p ajustada a FXH (eqs. 7-9)", "[V4]") {
    const v4::Config cf;
    const UnsatParams p = v4::Parametros(cf);
    for (double psi : {1., 3., 5., 7., 10., 15., 20.}) {
        REAL S, G;
        p.hydro.Retention(psi, S, G);
        const double C = 1. - std::log(1. + psi / 3000.) / std::log(1. + 1.e6 / 3000.);
        const double Sart = C / std::log(M_E + std::pow(psi / 5., 2.)); // theta_w / theta_s da eq. 7
        INFO("psi = " << psi << " kPa: S = " << S << ", eq. 7: " << Sart);
        CHECK(std::fabs(S - Sart) <= 0.002);
    }
    const double aev = 2.579;
    for (double psi : {3., 5., 7., 10., 12.5, 15., 20.}) {
        const auto st = p.hydro.Evaluate(-psi, 0.);
        const double fxh = KrFXH(psi, aev);
        INFO("psi = " << psi << " kPa: k_r = S^p = " << st.krw << ", FXH = " << fxh);
        CHECK(std::fabs(std::log10(st.krw / fxh)) <= 0.05);
    }
    const double psiP = v4::SuccaoPatamar(p.hydro, cf.qks);
    double lo = aev, hi = 100.;
    for (int i = 0; i < 60; ++i) { const double m = 0.5 * (lo + hi); (KrFXH(m, aev) > cf.qks ? lo : hi) = m; }
    WARN("patamar k_r = q/k_s: lei de potencia " << psiP << " kPa, FXH " << 0.5 * (lo + hi) << " kPa");
    CHECK(std::fabs(psiP / (0.5 * (lo + hi)) - 1.) <= 0.05);
}

// ---------------------------------------------------------------------------------------------------------------------
// Verificacao 2: poropressoes do regime permanente (q/k_s = 0,01) nas secoes X-X (x = 10 m) e Y-Y (x = 15 m) contra a
// Fig. 3 (dados/fig3_secao_*.csv, digitalizado, aproximado), malha de 1 m (o artigo usa 0,5 m; a de 0,5 m esta no teste
// [.lento]). Abaixo do lencol a poropressao e quase hidrostatica e depende so das cargas laterais e da recarga:
// tolerancia de 1,5 kPa (0,15 m de carga: a leitura, ~0,3 kPa, e a diferenca de malha). Acima do lencol o patamar de
// succao e o da funcao de permeabilidade: com a FXH, psi_aev desconhecido (de 0 a a = 5 kPa) move o patamar de 12,5 a
// 15,6 kPa (scripts/ajuste_kr.py), +-1,6 kPa em torno do nosso; com a leitura, tolerancia de 2 kPa. O lencol (p_w = 0)
// deve ficar a 0,15 m do do artigo. Confere tambem que o transiente chegou ao permanente e que o balanco de agua fecha.
TEST_CASE("V4: poropressoes permanentes nas secoes X-X e Y-Y (Fig. 3)", "[V4]") {
    const auto &d = CasoDeterministico("talude_liu_h10.msh");
    REQUIRE(d.perm.ok);
    INFO("chuva " << d.perm.chuva << ", esquerda " << d.perm.esq << ", direita " << d.perm.dir << ", face " << d.perm.face
                  << " m2/s");
    CHECK(std::fabs(d.perm.chuva + d.perm.esq + d.perm.dir + d.perm.face) <= 0.01 * d.perm.chuva);
    CHECK(std::fabs(d.perm.dV) <= 0.01 * d.perm.chuva);
    ComparaPerfil(d.xx, LeRef("fig3_secao_xx.csv"), 1.5, 2.0, "X-X");
    ComparaPerfil(d.yy, LeRef("fig3_secao_yy.csv"), 1.5, 2.0, "Y-Y");
}

// ---------------------------------------------------------------------------------------------------------------------
// Verificacao 3: FS deterministico. O artigo da FS = 1,525 pelo metodo de Bishop simplificado (SLOPE/W); o nosso e por
// reducao de resistencia em elementos finitos (SRM, Mohr-Coulomb associado), que difere do equilibrio limite: as
// comparacoes publicadas para taludes homogeneos c'-phi' (Griffiths e Lane 1999; Cheng, Lansivaara e Wei 2007) ficam
// dentro de ~5 %; a nossa succao no patamar e ~1,4 kPa menor (a FXH aproximada), o que baixa o FS em ~1 %. Faixa de +-6 %
// em torno de 1,525. O mecanismo deve ser uma ruptura de base (D_f entre 6 e 10 m; a superficie critica da Fig. 2 tem o
// ponto mais baixo em z = 2,9 m, D_f = 7,1 m) com a massa que desliza (|du| >= 0,5 max|du|) a 20 % da area entre o
// terreno e a superficie critica da Fig. 2 (44,5 m2, digitalizado; a area depende do limiar: README.md).
TEST_CASE("V4: FS deterministico por SRM contra Bishop (1,525)", "[V4]") {
    const auto &d = CasoDeterministico("talude_liu_h10.msh");
    WARN("FS (SRM) = " << d.fs.FS << " (Bishop 1,525), massa " << d.fs.massaArea << " m2, D_f = " << d.fs.Df << " m ("
                       << d.fs.modo << "), x de " << d.fs.massaXmin << " a " << d.fs.massaXmax << " m");
    CHECK(d.fs.FS >= 1.525 * 0.94);
    CHECK(d.fs.FS <= 1.525 * 1.06);
    CHECK(d.fs.modo == "base");
    CHECK(d.fs.Df <= 8.5);
    CHECK(std::fabs(d.fs.massaArea / 44.5 - 1.) <= 0.20);
}

// ---------------------------------------------------------------------------------------------------------------------
// Verificacao 4: os campos aleatorios de cada cenario (Tabela 3) e a sua estatistica (Tabela 1). s1: so k_s; s2: c' e
// phi'; s3: os tres. Comprimentos do nucleo exponencial: L = sqrt(pi) theta / 2 (26,59 m e 2,659 m), trocados com
// alfa = 90 graus. Em 400 realizacoes, a media e o desvio de ln c', ln phi' e ln k_mult num ponto devem ser os da
// lognormal (Tabela 1: COV 0,3, 0,2 e 0,6) dentro de 4 erros padrao (o desvio do campo truncado e sqrt(fracao da
// variancia de KL) vezes o alvo); com RHO = -0,5 (o valor do artigo) a correlacao entre ln c' e ln phi' no ponto deve
// ser -0,5 +- 0,15 (4 erros padrao de r com n = 400). A mesma SEMENTE e REAL refazem o mesmo campo.
TEST_CASE("V4: campos de KL por cenario (Tabelas 1 e 3)", "[V4]") {
    double lx, ly;
    v4::Config c0 = ConfigTeste("talude_liu_h10.msh");
    c0.Comprimentos(lx, ly);
    CHECK(lx == Approx(26.587).epsilon(1.e-3));
    CHECK(ly == Approx(2.6587).epsilon(1.e-3));
    c0.angulo = 90;
    c0.Comprimentos(lx, ly);
    CHECK(lx == Approx(2.6587).epsilon(1.e-3));
    CHECK(ly == Approx(26.587).epsilon(1.e-3));
    for (const char *caso : {"s1", "s2", "s3"}) {
        v4::Config cf = ConfigTeste("talude_liu_h10.msh");
        cf.caso = caso, cf.real = 3, cf.klM = 400;
        v4::Talude t(cf);
        t.AplicaCampos();
        const auto ec = v4::Estatistica(t, [&](size_t i) { return t.CPonto(i); });
        const auto ek = v4::Estatistica(t, [&](size_t i) { return t.KPonto(i); });
        INFO(caso << ": dp de c' " << ec.dp << ", dp de k_mult " << ek.dp);
        CHECK((ec.dp > 1.e-4 * ec.media) == (cf.caso != "s1")); // dp de um campo constante: so arredondamento
        CHECK((ek.dp > 1.e-4 * ek.media) == (cf.caso != "s2"));
        v4::Talude t2(cf); // mesma semente e realizacao
        t2.AplicaCampos();
        CHECK(t2.CPonto(100) == t.CPonto(100));
        CHECK(t2.KPonto(100) == t.KPonto(100));
    }
    // estatistica num ponto (o primeiro ponto de integracao com x > 15 m) em 400 realizacoes do cenario 3, RHO = -0,5
    v4::Config cf = ConfigTeste("talude_liu_h10.msh");
    cf.caso = "s3", cf.klM = 400, cf.rho = -0.5;
    const int N = 400;
    std::vector<double> lc(N), lp(N), lk(N);
    v4::Talude t(cf);
    size_t ip = 0;
    while (t.xs[ip] < 15.) ++ip;
    for (int r = 0; r < N; ++r) {
        t.cf.real = r;
        t.AplicaCampos();
        lc[r] = std::log(t.CPonto(ip)), lp[r] = std::log(t.PhiPonto(ip)), lk[r] = std::log(t.KPonto(ip));
    }
    // variancia do campo truncado no ponto (sum lambda_k f_k^2), para o desvio esperado
    const klsep::KL2D kl(0., 0., 30., 10., t.LxUsado, t.LyUsado, cf.klFrac, cf.klM);
    const double frac = kl.PreparaPontos({t.xs[ip]}, {t.ys[ip]}).var[0];
    auto media = [](const std::vector<double> &v) { double s = 0.; for (double x : v) s += x; return s / v.size(); };
    auto dp = [&](const std::vector<double> &v) { const double m = media(v); double s = 0.; for (double x : v) s += (x - m) * (x - m); return std::sqrt(s / (v.size() - 1)); };
    auto corr = [&](const std::vector<double> &a, const std::vector<double> &b) {
        const double ma = media(a), mb = media(b);
        double s = 0.;
        for (size_t i = 0; i < a.size(); ++i) s += (a[i] - ma) * (b[i] - mb);
        return s / (a.size() - 1) / (dp(a) * dp(b));
    };
    WARN("KL: " << 100. * t.klFracUsada << " % da variancia com 400 termos (" << 100. * frac << " % no ponto); dp ln c' " << dp(lc) << ", dp ln phi' "
                << dp(lp) << ", dp ln k " << dp(lk) << ", corr(ln c', ln phi') " << corr(lc, lp));
    for (auto [v, mu, cov] : {std::tuple{&lc, 10., 0.3}, std::tuple{&lp, 15., 0.2}, std::tuple{&lk, 1., 0.6}}) {
        const klsep::Lognormal L(mu, cov);
        const double s = L.sigma_ln * std::sqrt(frac), ep = s / std::sqrt(double(N));
        CHECK(std::fabs(media(*v) - L.mu_ln) <= 4. * ep);
        CHECK(std::fabs(dp(*v) / s - 1.) <= 4. / std::sqrt(2. * N));
    }
    CHECK(std::fabs(corr(lc, lp) + 0.5) <= 0.15);
    CHECK(std::fabs(corr(lc, lk)) <= 0.2);
}

// ---------------------------------------------------------------------------------------------------------------------
// Verificacao 5: o contrato do executavel com Monte Carlo (teste de fumaca). Uma realizacao do cenario 3 (alfa = 0) e uma
// do cenario 2 (alfa = 90, com o permanente lido de PW_PERMANENTE, pois k_s e constante), pelas variaveis de ambiente
// do contrato, na malha de 1 m e com 400 termos de KL: resultado.txt (real=, FS_min=, situacao=ok, seg=), resumo.csv
// (real e FS_min) e fs.csv (t_h e FS_SRM). O FS de uma realizacao deve ficar na faixa da Fig. 11 (minimo e maximo dos
// FS dos cenarios 2 e 3: ~1,0 a ~2,3).
TEST_CASE("V4: Monte Carlo, teste de fumaca do contrato do executavel", "[V4]") {
    const auto &d = CasoDeterministico("talude_liu_h10.msh");
    const auto base = std::filesystem::temp_directory_path() / "v4_fumaca";
    std::filesystem::remove_all(base);
    struct Caso { const char *caso, *angulo, *real, *pw; };
    for (const Caso &c : {Caso{"s3", "0", "1", ""}, Caso{"s2", "90", "2", "pw"}}) {
        setenv("CASO", c.caso, 1);
        setenv("ANGULO", c.angulo, 1);
        setenv("REAL", c.real, 1);
        setenv("CAMPOS", "1", 1);
        setenv("MALHA", "talude_liu_h10.msh", 1);
        setenv("KL_M", "400", 1);
        setenv("PW_PERMANENTE", *c.pw ? d.pw.c_str() : "", 1);
        const v4::Config cf = v4::Config::DoAmbiente();
        v4::Verboso() = 0;
        const auto dir = (base / (std::string(c.caso) + "_" + c.angulo)).string();
        const auto r = v4::RodaExecucao(cf, dir);
        CHECK(r.codigo == 0);
        CHECK(r.perm.lido == (*c.pw != 0));
        std::ifstream f(dir + "/resultado.txt");
        std::string l;
        REQUIRE(std::getline(f, l));
        INFO(c.caso << " alfa " << c.angulo << ": " << l);
        for (const char *k : {"real=", "FS_min=", "situacao=ok", "seg="}) CHECK(l.find(k) != std::string::npos);
        std::ifstream fr(dir + "/resumo.csv"), ff(dir + "/fs.csv");
        std::string cab, val;
        REQUIRE(std::getline(fr, cab));
        REQUIRE(std::getline(fr, val));
        CHECK(cab.rfind("real,", 0) == 0);
        CHECK(cab.find(",FS_min,") != std::string::npos);
        REQUIRE(std::getline(ff, cab));
        CHECK(cab.rfind("t_h,FS_SRM", 0) == 0);
        CHECK(r.fs.FS >= 1.0);
        CHECK(r.fs.FS <= 2.3);
    }
    for (const char *k : {"CASO", "ANGULO", "REAL", "CAMPOS", "MALHA", "KL_M", "PW_PERMANENTE"}) unsetenv(k);
}

// ---------------------------------------------------------------------------------------------------------------------
// [.lento] Verificacoes 2 e 3 na malha do artigo (0,5 m): mesmas tolerancias e faixas (ver acima).
TEST_CASE("V4: secoes e FS na malha de 0,5 m do artigo", "[V4][.lento]") {
    const auto &d = CasoDeterministico("talude_liu_h05.msh");
    REQUIRE(d.perm.ok);
    ComparaPerfil(d.xx, LeRef("fig3_secao_xx.csv"), 1.5, 2.0, "X-X (h = 0,5 m)");
    ComparaPerfil(d.yy, LeRef("fig3_secao_yy.csv"), 1.5, 2.0, "Y-Y (h = 0,5 m)");
    WARN("h = 0,5 m: FS (SRM) = " << d.fs.FS << ", massa " << d.fs.massaArea << " m2, D_f = " << d.fs.Df << " m");
    CHECK(d.fs.FS >= 1.525 * 0.94);
    CHECK(d.fs.FS <= 1.525 * 1.06);
    CHECK(d.fs.modo == "base");
    CHECK(std::fabs(d.fs.massaArea / 44.5 - 1.) <= 0.20);
}

// [.lento] Monte Carlo curto (N = V4_N_MC realizacoes, padrao 16, malha de 1 m) dos cenarios 1 (k_s) e 2 (c' e phi', com
// RHO = -0,5 como no artigo e o permanente do deterministico) a alfa = 0, contra a Fig. 11 e a Tabela 4. Media: a do
// artigo (1,52 e 1,49, calculadas da Tabela 4 e da Fig. 11) corrigida pelo vies deterministico SRM/Bishop (FS_det/1,525),
// dentro de 3 erros padrao + 0,02. Desvio: o nucleo exponencial (no lugar do gaussiano, com a mesma escala de flutuacao)
// e o erro amostral (~18 % com N = 16) deixam uma faixa de fator 2 (0,5 a 2 vezes 0,041 e 0,149). Os P_f do artigo
// (0 e 4e-5) nao podem ser resolvidos com N pequeno: so o beta lognormal (eq. 13) e informado.
TEST_CASE("V4: Monte Carlo curto dos cenarios 1 e 2 (Fig. 11)", "[V4][.lento]") {
    const auto &d = CasoDeterministico("talude_liu_h10.msh");
    const int N = v4::EnvI("V4_N_MC", 16);
    const auto base = std::filesystem::temp_directory_path() / "v4_mc_curto";
    struct Cen { const char *caso; double rho, media, dp; };
    for (const Cen &c : {Cen{"s1", 0., 1.52, 0.041}, Cen{"s2", -0.5, 1.49, 0.149}}) {
        std::vector<double> fs;
        for (int r = 1; r <= N; ++r) {
            v4::Config cf = ConfigTeste("talude_liu_h10.msh");
            cf.caso = c.caso, cf.real = r, cf.rho = c.rho, cf.klM = 1000;
            if (std::string(c.caso) == "s2") cf.pwPermanente = d.pw;
            const auto res = v4::RodaExecucao(cf, (base / (std::string(c.caso) + "_r" + std::to_string(r))).string());
            REQUIRE(res.codigo == 0);
            fs.push_back(res.fs.FS);
        }
        double m = 0., s = 0.;
        for (double v : fs) m += v;
        m /= N;
        for (double v : fs) s += (v - m) * (v - m);
        s = std::sqrt(s / (N - 1));
        const double cov = s / m, beta = std::log(m / std::sqrt(1. + cov * cov)) / std::sqrt(std::log(1. + cov * cov));
        const double alvo = c.media * d.fs.FS / 1.525;
        WARN(c.caso << " (rho = " << c.rho << "), N = " << N << ": media FS " << m << " (artigo " << c.media
                    << ", corrigida pelo vies deterministico " << alvo << "), dp " << s << " (artigo " << c.dp
                    << "), beta lognormal " << beta);
        CHECK(std::fabs(m - alvo) <= 3. * s / std::sqrt(double(N)) + 0.02);
        CHECK(s >= 0.5 * c.dp);
        CHECK(s <= 2.0 * c.dp);
    }
}
