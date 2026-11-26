// Kang2020: a reproducao do exemplo de Kang, Cho, Kim e Go (2020), Water 12(3):812, escrita passo a passo.
// O exemplo (Secoes 2.3 e 3.2 do artigo): um talude infinito com 2 m de solo (medidos na vertical) sobre uma base
// impermeavel, sob a chuva horaria de 26-27/07/2011 no morro Umyeon (Seul), durante 17 h (das 16:00 as 09:00), com
// Pw inicial de -20 kPa constante com a profundidade e Pa = 0, nas 5 zonas da Tabela 1. O artigo compara o fluxo so da
// agua (ar passivo, "single-phase", SP) com o fluxo bifasico agua-ar (ar ativo, "two-phase", TP), no FLAC 2D acoplado,
// e calcula o FS do talude infinito (Eq. 14), o minimo ao longo da profundidade (Figs. 7, 8 e 9).
//
// Uso: Kang2020 [passo=1|2|3] [zona=1..5] [modelo=sp|tp] [ny=40] [dt=0.05] [saida=kang] [vgtab=auto|0|<eps>]
//               [kramin=1e-4] [invP0=<1/kPa>] [malha=<arquivo.csv>] [geometria=coluna|faixa|normal] [beta=<graus>]
//               [chuva_h=0|1]
//   malha: so monta a coluna, grava a malha geometrica (os triangulos do solo e as linhas das fronteiras, com o id de
//   cada uma) e sai; figura_malha.py desenha
// Passo 6: a malha do artigo (geometria=faixa). A faixa do FLAC da Fig. 5: um paralelogramo com 10 m na horizontal,
//   laterais verticais de 2 m e o angulo beta (o padrao: o angulo medio da zona; o artigo desenha 31 graus), na grade
//   de 40 x 20 zonas (0,25 m na horizontal, 0,1 m na vertical), cada zona em 4 triangulos em volta do centro (os nos
//   do FLAC mais os centros). Gravidade vertical: a agua escorre ao longo da faixa e se acumula embaixo, entao o fluxo
//   depende do angulo (uma rodada por angulo). Contornos da Fig. 5: topo com q_w = I(t) (normal a superficie
//   inclinada, como no FLAC) e p_a = 0 (TP); laterais u_x = 0; base u = 0; sem fluxo nas laterais e na base. O perfil
//   do FS e a vertical do meio (x = 5 m), e o FS da Eq. 14 sai no angulo da faixa (FS_beta)
// Passo 7: a coluna normal ao talude (geometria=normal), o talude infinito 1D de verdade. O eixo da coluna e a normal
//   a superficie; a gravidade que importa ao longo dele e g cos(beta) (SetGravityDirection(0, -cos(beta)) no material:
//   escala o peso do solo e os pesos da agua e do ar nos termos de Darcy). A componente paralela, g sin(beta), fica de
//   fora: no talude infinito o fluxo paralelo e uniforme (sem divergencia) e nao muda o perfil de pressoes. A coluna
//   tem a espessura normal H cos(beta) (os 2 m verticais do artigo), a chuva entra como I cos(beta) por metro de
//   superficie (a chuva vertical, medida por metro horizontal) e o perfil volta para a profundidade vertical
//   z = d / cos(beta) na Eq. 14. Saturada sobre a base impermeavel, ela da p_w = gamma_w d cos(beta) =
//   gamma_w z cos^2(beta), o escoamento paralelo, como o meio da faixa. chuva_h=1: na faixa, a chuva por metro horizontal
//   (q = I cos(beta) normal a superficie inclinada) no lugar de I por metro inclinado
//
// Passo 1: os dados do artigo (kang_dados.h) e duas conferencias, sem transiente:
//   1. as relacoes hidraulicas da biblioteca (TPZUnsaturatedHydraulics: van Genuchten, Mualem da agua e do ar, as
//      mesmas Eqs. 1, 4 e 5 do artigo) contra as curvas da Fig. 4;
//   2. o FS inicial pela Eq. 14 (talude_infinito.h), com Pw = -20 kPa e Pa = 0 em todo o perfil, contra o t = 0 da
//      Fig. 9, nos tres angulos de cada zona.
// Passo 2 (este): a coluna no neopz, ainda sem chuva:
//   1. a malha: a faixa [0, w] x [0, H] (MalhaColuna, a mesma do SlopeRainColumn), triangulos em volta de um no central
//      em cada celula, simetrica no espelho x -> w - x;
//   2. o material u - p_w - p_a da biblioteca (TPZMatUnsaturatedUPwPa) com esqueleto elastico e as equacoes de massa
//      acopladas a deformacao (o FLAC do artigo e acoplado), a curva de van Genuchten - Mualem da zona, chi = S_w;
//   3. os contornos: base engastada, laterais presas so na horizontal, e tudo impermeavel menos o topo (sem fluxo e
//      a condicao natural); no topo a chuva, e no modelo TP a pressao de ar atmosferica (Pa = 0, a "atmospheric
//      boundary" do artigo; o resto e estanque ao ar);
//   4. o estado inicial: succao uniforme de 20 kPa (Pw = -20 kPa, Pa = 0) e o estado geostatico (RainTransient);
//   5. 17 h sem chuva, para conferir o modelo antes de molhar: a agua guardada tem de se conservar (as fronteiras sao
//      fechadas), a faixa tem de ficar simetrica (o problema e 1D) e a agua so se redistribui, descendo para a base.
// Passo 3 (este): a chuva da Fig. 2 hora a hora (SetIntensity e AdvanceTo ate o fim de cada hora). No topo, o
//   RainTransient aplica o fluxo da chuva enquanto a agua entra e troca o no por p_w = 0 (a face de percolacao: o
//   excesso escorre) quando ela empocaria. A cada hora, o perfil (p_w, p_a, S_w) no eixo da faixa e o FS da Eq. 14 nos
//   tres angulos da zona; os CSV <saida>_perfis e <saida>_fs vao para a comparacao com as Figs. 7, 8 e 9
//   (compara_kang.py). Rodar com modelo=sp compara com a Fig. 7 e as curvas SP da Fig. 9.
//   vgtab: no solo do norte (n_VG = 1,54 < 2) a derivada de k_r de Mualem explode na saturacao e o Newton cicla quando
//   a superficie satura (a zona 2 parava em 11,8 h). O remedio e o do SlopeRainColumn (Cho 2014): k_r tabelado
//   (MualemTabelado.h) ate S_e = 1 - eps e fechado suavemente ate k_r = 1 na saturacao. auto (o padrao): so com n < 2
// Passo 4: o mesmo Passo 3 com modelo=tp (o bifasico: p_a livre, com p_a = 0 so nos nos do topo e o resto estanque ao
//   ar), comparado com a Fig. 8 e as curvas TP da Fig. 9. kramin: o piso da permeabilidade relativa do ar (o
//   fKraMin da biblioteca, 1e-4 por padrao). Com o piso, o ar ainda atravessa a zona saturada (k_a = 1e-4 x 56 x k_s,
//   ~0,6 % do k_s), o que pode aliviar o p_a preso embaixo da lamina d'agua; o artigo (Eq. 5) vai a zero

#include "kang_dados.h"
#include "talude_infinito.h"

// as classes do SlopeRandomFieldRefinement (o modelo, o material com hp e o transiente da chuva)
#include "HPMaterial.h"
#include "MualemTabelado.h"
#include "RainTransient.h"
#include "UnsatModel.h"

// o material nao saturado u - p_w - p_a da biblioteca, com Mohr-Coulomb
#include "Plasticity/TPZElasticResponse.h"
#include "Plasticity/TPZMatUnsaturatedUPwPa.h"
#include "Plasticity/TPZPlasticStepVoigt.h"
#include "Plasticity/TPZYCMohrCoulombPV2.h"

// a malha geometrica
#include "TPZGeoLinear.h"
#include "pzgeotriangle.h"
#include "tpzgeoelrefpattern.h"

#include <array>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <map>
#include <memory>
#include <set>
#include <string>
#include <vector>

using namespace kang;

using TMC = TPZPlasticStepVoigt<TPZYCMohrCoulombPV2, TPZElasticResponse>;
using TUnsat = TPZMatUnsaturatedUPwPa<TMC, TPZUnsaturatedMem>;
using TUnsatHP = HPMaterial<TUnsat>;

// os ids da malha: o solo e as quatro fronteiras (os do SlopeRainColumn)
enum EIds { matSolo = 1, base = -1, direita = -2, esquerda = -5, topo = -6 };

// ---- as opcoes da linha de comando (chave=valor)
struct Opcoes {
    int passo = 2, zona = 1, ny = 40;
    bool tp = false;   // modelo: false = SP (so a agua, ar passivo), true = TP (bifasico, ar ativo)
    REAL dt = 0.05;    // passo de tempo (h)
    std::string saida = "kang"; // prefixo dos CSV (Passo 3)
    REAL vgtab = -1.;           // eps do k_r tabelado: < 0 = auto (1e-3 se n_VG < 2), 0 = Mualem direto
    REAL invP0 = 0.;            // > 0: troca o 1/P0 (1/kPa) da curva da zona (teste da curva do norte da Fig. 4a)
    REAL kramin = 1.e-4;        // piso de k_ra (Passo 4)
    std::string malha;          // grava a malha neste CSV e sai
    bool faixa = false;         // geometria: false = a coluna vertical (Passos 2 a 4), true = a faixa do FLAC (Passo 6)
    bool normal = false;        // geometria=normal: a coluna normal ao talude (Passo 7)
    bool chuvaH = false;        // chuva_h=1: na faixa, a chuva por metro horizontal (Passo 7)
    REAL beta = -1.;            // angulo da faixa (graus); < 0: o angulo medio da zona
    int nx = 40;                // zonas da faixa ao longo dos 10 m (a espessura usa ny; o artigo: 40 x 20)
};

static Opcoes LeOpcoes(int argc, char **argv) {
    Opcoes op;
    bool nyDado = false;
    for (int i = 1; i < argc; i++) {
        const std::string a = argv[i];
        const size_t eq = a.find('=');
        const std::string k = a.substr(0, eq), v = eq == std::string::npos ? "" : a.substr(eq + 1);
        if (k == "passo") op.passo = std::stoi(v);
        else if (k == "zona") op.zona = std::stoi(v);
        else if (k == "ny") op.ny = std::stoi(v), nyDado = true;
        else if (k == "dt") op.dt = std::stod(v);
        else if (k == "saida") op.saida = v;
        else if (k == "vgtab") op.vgtab = v == "auto" ? -1. : std::stod(v);
        else if (k == "invP0") op.invP0 = std::stod(v);
        else if (k == "kramin") op.kramin = std::stod(v);
        else if (k == "malha") op.malha = v;
        else if (k == "geometria" && (v == "coluna" || v == "faixa" || v == "normal"))
            op.faixa = v == "faixa", op.normal = v == "normal";
        else if (k == "chuva_h") op.chuvaH = std::stoi(v) != 0;
        else if (k == "beta") op.beta = std::stod(v);
        else if (k == "nx") op.nx = std::stoi(v);
        else if (k == "modelo" && (v == "sp" || v == "tp")) op.tp = v == "tp";
        else {
            std::fprintf(stderr, "opcao desconhecida: %s\nuso: Kang2020 [passo=1|2|3] [zona=1..5] [modelo=sp|tp] "
                                 "[ny=40] [dt=0.05] [saida=kang] [vgtab=auto|0|<eps>]\n", argv[i]);
            std::exit(1);
        }
    }
    if (op.zona < 1 || op.zona > 5) std::exit((std::fprintf(stderr, "zona: 1 a 5\n"), 1));
    if (op.faixa && !nyDado) op.ny = 20; // a grade do artigo: 40 x 20 zonas
    return op;
}

// ======================================================================================================== Passo 1
// As curvas hidraulicas de um lado do morro numa tabela: a succao, o teor de umidade theta = n S_w, S_w, S_e e as
// permeabilidades relativas da agua e do ar
static void TabelaHidraulica(const Lado &l) {
    const TPZUnsaturatedHydraulics h = l.Hidraulica();
    std::printf("\n  lado %s: 1/P0 = %.3f 1/kPa, m = %.2f (n_VG = %.3f), theta_s = %.3f, theta_r = %.3f (S_rw = %.3f)\n",
                l.nome, l.invP0, l.m, h.fBVG, l.thetaS, l.thetaR, l.Srw());
    std::printf("  %9s %7s %7s %7s %10s %10s\n", "psi (kPa)", "theta", "S_w", "S_e", "k_rw", "k_ra");
    for (REAL psi : {0.1, 1., 5., 10., 20., 50., 100., 1000., 1.e5}) {
        REAL S, Gamma, Se, dSe, krw, dkrw, kra, dkra;
        h.Retention(psi, S, Gamma); // S_w(psi) e a derivada (nao usada aqui)
        h.Effective(S, Se, dSe);    // S_e = (S_w - S_rw)/(1 - S_rw)
        h.WaterPermeability(S, krw, dkrw);
        h.AirPermeability(S, kra, dkra);
        std::printf("  %9g %7.3f %7.3f %7.3f %10.3e %10.3e\n", psi, l.thetaS * S, S, Se, krw, kra);
    }
}

static void Passo1() {
    std::printf("Kang et al. (2020), Passo 1: os dados e as conferencias iniciais\n");
    std::printf("\n1. Relacoes hidraulicas (compare com a Fig. 4)\n");
    TabelaHidraulica(Norte());
    std::printf("  LIDO da Fig. 4 (norte): theta = 0,44 / 0,40 / 0,30 / 0,235 a 10 / 20 / 100 / 1000 kPa; "
                "k_r = 0,13 / 0,05 / 4e-4 a 10 / 20 / 100 kPa\n");
    TabelaHidraulica(Sul());
    std::printf("  LIDO da Fig. 4 (sul): theta = 0,49 / 0,43 / 0,21 / 0,185 a 10 / 20 / 100 / 1000 kPa; "
                "k_r = 0,8 / 0,25 / 1e-4 a 10 / 20 / 100 kPa\n");

    std::printf("\n2. FS inicial pela Eq. 14 (Pw = %g kPa, Pa = 0), o minimo ao longo dos %g m; compare com t = 0 da "
                "Fig. 9\n", pw0, H);
    std::printf("  %4s %6s %6s %8s %8s | %26s | %26s | %20s\n", "zona", "n", "S_w0", "gEq15", "gamma_t",
                "FS (W da Eq. 15) max/med/min", "FS (W = gamma_t z)", "Fig. 9 (lido)");
    for (const Zona &zo : Zonas()) {
        const Lado &l = zo.L();
        const TPZUnsaturatedHydraulics h = l.Hidraulica();
        const REAL n = l.thetaS;
        REAL Sw0, Gamma;
        h.Retention(-pw0, Sw0, Gamma); // a succao inicial e Pa - Pw = 20 kPa
        std::vector<PontoPerfil> perfil; // o perfil uniforme de 0 a H, de 1 em 1 cm
        for (int i = 0; i <= 200; i++) perfil.push_back({H * i / 200., Sw0, pw0, 0.});
        const REAL gEq15 = zo.gammaD + n * Sw0 * gammaW; // o peso especifico da Eq. (15) no estado inicial
        REAL fsA[3], fsB[3];
        for (int k = 0; k < 3; k++) {
            const REAL b = zo.beta[k] * M_PI / 180., phi = zo.phi * M_PI / 180.;
            fsA[k] = FSInfinito(perfil, b, phi, zo.c, zo.gammaD, n, gammaW);
            fsB[k] = FSInfinito(perfil, b, phi, zo.c, zo.gammaT, 0., gammaW); // W = gamma_t z (gamma_d -> gamma_t, n = 0)
        }
        std::printf("  %4d %6.3f %6.3f %8.2f %8.2f | %8.3f %8.3f %8.3f | %8.3f %8.3f %8.3f | %6.2f %6.2f %6.2f\n",
                    zo.id, n, Sw0, gEq15, zo.gammaT, fsA[0], fsA[1], fsA[2], fsB[0], fsB[1], fsB[2], zo.fs0Fig9[0],
                    zo.fs0Fig9[1], zo.fs0Fig9[2]);
    }
    std::printf("\n3. A chuva da Fig. 2 (digitalizada), das 16:00 de 26/07 as 09:00 de 27/07/2011 (mm/h):\n  ");
    REAL total = 0.;
    for (size_t i = 0; i < Chuva().size(); i++) std::printf("%s%.1f", i ? ", " : "", Chuva()[i]), total += Chuva()[i];
    std::printf("\n  total %.1f mm em %g h\n", total, duracao);
}

// ======================================================================================================== Passo 2
// 1. A malha da coluna: a faixa [0, w] x [0, H] com ny celulas de altura H/ny. Nos dos cantos: 2j (x = 0) e 2j + 1
// (x = w) na cota y_j = H j / ny; centros: 2 (ny + 1) + j. Cada celula vira 4 triangulos em volta do centro (assim a
// malha e simetrica e o problema fica 1D de verdade). As fronteiras ganham elementos lineares com os ids de EIds: e por
// eles que a analise acha os apoios de u e o RainTransient acha os nos do topo, que recebem a chuva
static TPZGeoMesh *MalhaColuna(REAL altura, REAL w, int ny) {
    auto *g = new TPZGeoMesh();
    g->SetDimension(2);
    const int64_t nCantos = 2 * int64_t(ny + 1);
    g->NodeVec().Resize(nCantos + ny);
    auto poe = [g](int64_t id, REAL x, REAL y) {
        TPZManVector<REAL, 3> c = {x, y, 0.};
        g->NodeVec()[id] = TPZGeoNode(id, c, *g);
    };
    auto Y = [altura, ny](int j) { return j == ny ? altura : altura * j / ny; };
    for (int j = 0; j <= ny; j++) poe(2 * j, 0., Y(j)), poe(2 * j + 1, w, Y(j));
    for (int j = 0; j < ny; j++) poe(nCantos + j, 0.5 * w, 0.5 * (Y(j) + Y(j + 1)));
    for (int j = 0; j < ny; j++) {
        const int64_t bl = 2 * j, br = 2 * j + 1, tl = 2 * j + 2, tr = 2 * j + 3, c = nCantos + j;
        for (const auto &t : {std::array<int64_t, 3>{bl, br, c}, std::array<int64_t, 3>{br, tr, c},
                              std::array<int64_t, 3>{tr, tl, c}, std::array<int64_t, 3>{tl, bl, c}}) {
            TPZManVector<int64_t, 3> v = {t[0], t[1], t[2]};
            new TPZGeoElRefPattern<pzgeom::TPZGeoTriangle>(v, matSolo, *g); // anti-horarios
        }
    }
    auto linha = [g](int64_t a, int64_t b, int id) {
        TPZManVector<int64_t, 2> v = {a, b};
        new TPZGeoElRefPattern<pzgeom::TPZGeoLinear>(v, id, *g);
    };
    linha(0, 1, base);
    for (int j = 0; j < ny; j++) linha(2 * j + 1, 2 * j + 3, direita);
    linha(2 * int64_t(ny) + 1, 2 * int64_t(ny), topo);
    for (int j = ny; j > 0; j--) linha(2 * j, 2 * j - 2, esquerda);
    g->BuildConnectivity();
    return g;
}

// Passo 6: a faixa do FLAC (Fig. 5 do artigo). Vertices (i, j), i = 0..nx ao longo de x e j = 0..ny na espessura:
// x_i = L i / nx, y = x_i tan(beta) + H j / ny (a base e a reta y = x tan(beta), a superficie H acima dela), no id
// i (ny + 1) + j; os centros das zonas depois deles. Cada zona (um paralelogramo de laterais verticais) em 4 triangulos
// anti-horarios em volta do centro; as fronteiras com os mesmos ids da coluna
static TPZGeoMesh *MalhaFaixa(REAL L, REAL altura, REAL beta, int nx, int ny) {
    auto *g = new TPZGeoMesh();
    g->SetDimension(2);
    const int64_t nv = int64_t(nx + 1) * (ny + 1);
    g->NodeVec().Resize(nv + int64_t(nx) * ny);
    auto id = [ny](int i, int j) { return int64_t(i) * (ny + 1) + j; };
    auto poe = [g](int64_t k, REAL x, REAL y) {
        TPZManVector<REAL, 3> c = {x, y, 0.};
        g->NodeVec()[k] = TPZGeoNode(k, c, *g);
    };
    const REAL tb = std::tan(beta);
    for (int i = 0; i <= nx; i++)
        for (int j = 0; j <= ny; j++) {
            const REAL x = L * i / nx;
            poe(id(i, j), x, x * tb + altura * j / ny);
        }
    for (int i = 0; i < nx; i++)
        for (int j = 0; j < ny; j++) {
            const REAL x = L * (i + 0.5) / nx;
            poe(nv + int64_t(i) * ny + j, x, x * tb + altura * (j + 0.5) / ny);
            const int64_t bl = id(i, j), br = id(i + 1, j), tr = id(i + 1, j + 1), tl = id(i, j + 1);
            const int64_t c = nv + int64_t(i) * ny + j;
            for (const auto &t : {std::array<int64_t, 3>{bl, br, c}, std::array<int64_t, 3>{br, tr, c},
                                  std::array<int64_t, 3>{tr, tl, c}, std::array<int64_t, 3>{tl, bl, c}}) {
                TPZManVector<int64_t, 3> v = {t[0], t[1], t[2]};
                new TPZGeoElRefPattern<pzgeom::TPZGeoTriangle>(v, matSolo, *g);
            }
        }
    auto linha = [g](int64_t a, int64_t b, int k) {
        TPZManVector<int64_t, 2> v = {a, b};
        new TPZGeoElRefPattern<pzgeom::TPZGeoLinear>(v, k, *g);
    };
    for (int i = 0; i < nx; i++) linha(id(i, 0), id(i + 1, 0), base);         // a base, da esquerda para a direita
    for (int j = 0; j < ny; j++) linha(id(nx, j), id(nx, j + 1), direita);     // a lateral direita, subindo
    for (int i = nx; i > 0; i--) linha(id(i, ny), id(i - 1, ny), topo);        // o topo, da direita para a esquerda
    for (int j = ny; j > 0; j--) linha(id(0, j), id(0, j - 1), esquerda);      // a lateral esquerda, descendo
    g->BuildConnectivity();
    return g;
}

// 2. O material da zona: u - p_w - p_a com esqueleto elastico. O Mohr-Coulomb vai com uma coesao enorme, para a
// projecao ser sempre a elastica (o transiente so da as pressoes e a deformacao; o FS vem da Eq. 14, no pos-processamento)
static TUnsatHP *MaterialKang(const Zona &z, REAL vgtab, REAL invP0, REAL kramin) {
    const REAL phi = z.phi * M_PI / 180.;
    TPZElasticResponse ER;
    ER.SetEngineeringData(E, nu);
    auto *mat = new TUnsatHP(matSolo);
    mat->SetPlasticModel(TMC(TPZYCMohrCoulombPV2(phi, phi, 1.e9, ER), ER));
    TPZUnsaturatedHydraulics h = z.L().Hidraulica(); // van Genuchten - Mualem da agua e do ar (Passo 1)
    if (invP0 > 0.) h.fAVG = invP0;                   // o teste da curva (opcao invP0)
    h.fKraMin = kramin;                               // o piso de k_ra (opcao kramin)
    const REAL eps = vgtab < 0. ? (h.fBVG < 2. ? 1.e-3 : 0.) : vgtab;
    if (eps > 0.) // k_r da agua tabelado ate S_e = 1 - eps (Passo 3, vgtab)
        UsaTabela(h, MualemTabelado(h.fBVG, h.fSwr, h.fAVG, eps, EFechoVG::C1));
    mat->SetHydraulics(h);
    mat->SetBishop(TUnsat::EChiSaturation);       // chi = S_w (Eq. 12 do artigo)
    mat->SetPorosity(z.L().thetaS);               // n = theta_s (Fig. 4a)
    mat->SetBiot(1., 0.);
    mat->SetWaterBulkModulus(Kw);                 // Tabela 2
    mat->SetConductivity(z.ks, gammaW, muRatio);  // k_s da zona; o ar com mu_w / mu_a = 56 (Tabela 2)
    mat->SetAir(rhoA * 9.81e-3, 101.325);         // peso especifico do ar (kN/m3) e a pressao atmosferica (kPa)
    mat->SetUnitWeight(z.gammaT);                 // o peso no equilibrio (o FS usa a Eq. 15)
    mat->SetLumpedStorage(true);                  // armazenamento concentrado nos vertices (estavel com P1)
    mat->SetHydroMechanicalCoupling(true);        // os termos S alfa d(eps_v)/dt nas equacoes de massa (acoplado)
    return mat;
}

// 4. Os dados do transiente: succao inicial uniforme de 20 kPa (o nivel d'agua bem fundo e o teto de succao de 20 kPa:
// p_w = -min(gamma_w x altura acima do nivel, 20 kPa) = -20 kPa em todo ponto), sem chuva por enquanto
static RainData DadosKang(const Opcoes &op) {
    RainData r;
    r.intensity = 0.;      // mm/h; o driver troca hora a hora (SetIntensity)
    r.duration = 1.e9;     // a chuva nao acaba sozinha: quem decide e o driver
    r.dt = op.dt;          // h
    r.dtMin = op.dt / 64.; // um passo que nao converge e refeito em metades, ate dt/64
    r.lineSearch = 1;      // o Newton do SlopeRainColumn (newton=1e-6,80,1)
    r.rainVertical = op.chuvaH; // coluna: o topo e horizontal (as duas convencoes dao o mesmo); faixa: chuva_h
    r.xToe = 0., r.Hw = 0., r.betaW = 0., r.yToe = -1.e6; // o nivel d'agua 1000 km abaixo: so o teto de succao vale
    r.suctionMax = -pw0;
    return r;
}

// O perfil no eixo esquerdo da faixa (x = 0), de cima para baixo: os nos 2j com z = H - y_j
struct Perfil {
    std::vector<int64_t> nos;
    std::vector<REAL> z;
};
static Perfil FazPerfil(int ny) {
    Perfil P;
    for (int j = ny; j >= 0; j--) P.nos.push_back(2 * j), P.z.push_back(H - H * j / ny);
    return P;
}
// Passo 6: na faixa, a vertical do meio (os vertices i = nx/2), de cima para baixo
static Perfil FazPerfilFaixa(int nx, int ny) {
    Perfil P;
    const int i = nx / 2;
    for (int j = ny; j >= 0; j--) P.nos.push_back(int64_t(i) * (ny + 1) + j), P.z.push_back(H - H * j / ny);
    return P;
}

// A agua guardada na faixa (m3 por metro de espessura): o mesmo armazenamento concentrado do material, a area de cada
// triangulo dividida pelos tres vertices, n S_w(p_a - p_w) em cada um
static REAL AguaGuardada(TPZGeoMesh *g, const std::vector<REAL> &pw, const std::vector<REAL> &pa,
                         const TPZUnsaturatedHydraulics &h, REAL n) {
    REAL soma = 0.;
    for (TPZGeoEl *gel : g->ElementVec()) {
        if (!gel || gel->Dimension() != 2) continue;
        std::array<std::array<REAL, 2>, 3> x;
        for (int i = 0; i < 3; i++) x[i] = {gel->NodePtr(i)->Coord(0), gel->NodePtr(i)->Coord(1)};
        const REAL area = 0.5 * std::fabs((x[1][0] - x[0][0]) * (x[2][1] - x[0][1]) -
                                          (x[2][0] - x[0][0]) * (x[1][1] - x[0][1]));
        for (int i = 0; i < 3; i++) {
            const int64_t no = gel->NodeIndex(i);
            REAL S, Gm;
            h.Retention(pa[no] - pw[no], S, Gm);
            soma += area / 3. * n * S;
        }
    }
    return soma;
}

// A coluna montada (os passos 1 a 4 acima): o modelo (malhas e material), a analise e o transiente da chuva, nesta
// ordem, porque cada um aponta para o anterior (e sao destruidos na ordem inversa)
struct Coluna {
    const Zona &zo;
    const Opcoes &op;
    UnsatModel<TUnsatHP> m;
    std::unique_ptr<TPZUnsatAnalysis> an;
    std::unique_ptr<RainTransient> chuva;
    Perfil P;

    REAL beta = 0.; // o angulo da faixa (rad; 0 na coluna)

    Coluna(const Zona &z, const Opcoes &o) : zo(z), op(o) {
        if (op.faixa || op.normal) beta = (op.beta < 0. ? zo.beta[1] : op.beta) * M_PI / 180.;
        if (op.faixa)
            m.SetGeoMesh(MalhaFaixa(Lfaixa, H, beta, op.nx, op.ny));
        else if (op.normal) // a coluna normal: a espessura normal H cos(beta), celulas quadradas
            m.SetGeoMesh(MalhaColuna(H * std::cos(beta), H * std::cos(beta) / op.ny, op.ny));
        else
            m.SetGeoMesh(MalhaColuna(H, H / op.ny, op.ny)); // celulas quadradas
        // os apoios de u (base engastada, laterais so em x); sem fluxo e a condicao natural (nada a declarar)
        m.AddDirichletBC(base, TUnsat::EDirichletU, 0, {0, 1});
        m.AddDirichletBC(esquerda, TUnsat::EDirichletU, 0, {0});
        m.AddDirichletBC(direita, TUnsat::EDirichletU, 0, {0});
        TUnsatHP *mat = MaterialKang(zo, op.vgtab, op.invP0, op.kramin);
        if (op.normal) // Passo 7: so a componente da gravidade normal ao talude, ao longo da coluna
            mat->SetGravityDirection(TPZManVector<REAL, 3>{0., -std::cos(beta), 0.});
        m.SetMaterial(mat);
        m.CreateCompMesh({{2, 2}, {1, 1}, {1, 1}}); // u P2, p_w P1, p_a P1 (Taylor-Hood)
        an = std::make_unique<TPZUnsatAnalysis>(
            m.CompMesh(), std::map<int, std::array<bool, 2>>{{base, {true, true}}, {direita, {true, false}},
                                                             {esquerda, {true, false}}});
        an->SetTolerances(1.e-7, 1.e-6, 80, 25);
        // o transiente: chuva no topo e, no TP, p_a = 0 nos nos do topo
        chuva = std::make_unique<RainTransient>(*an, DadosKang(op), gammaW, std::set<int>{}, std::set<int>{topo},
                                                std::map<int, REAL>{},
                                                op.tp ? ERainModel::ECoupledActiveAir : ERainModel::ECoupledPassiveAir);
        P = op.faixa ? FazPerfilFaixa(op.nx, op.ny) : FazPerfil(op.ny);
        if (op.normal) // a profundidade normal: a malha tem a altura H cos(beta)
            for (REAL &d : P.z) d *= std::cos(beta);
    }
    static constexpr REAL Lfaixa = 10.; // comprimento horizontal da faixa (m)
    // a largura horizontal que recebe a chuva (a agua infiltrada em mm e por metro horizontal); na coluna normal, a
    // largura w ao longo do talude projeta w cos(beta) na horizontal
    REAL Largura() const {
        if (op.faixa) return Lfaixa;
        if (op.normal) return H * std::cos(beta) / op.ny * std::cos(beta);
        return H / op.ny;
    }
    // a profundidade vertical de um ponto do perfil (na coluna normal o perfil e medido na normal)
    REAL Vertical(REAL d) const { return op.normal ? d / std::cos(beta) : d; }
    const TPZUnsaturatedHydraulics &Hidraulica() { return m.Material()->Hydraulics(); }
    // o perfil do eixo esquerdo agora: z, S_w, p_w, p_a (o PontoPerfil da Eq. 14)
    std::vector<PontoPerfil> PerfilAgora() {
        const std::vector<REAL> pw = chuva->NodalWaterPressure(), pa = chuva->NodalAirPressure();
        std::vector<PontoPerfil> v;
        for (size_t i = 0; i < P.nos.size(); i++) {
            const int64_t no = P.nos[i];
            REAL S, Gm;
            Hidraulica().Retention(pa[no] - pw[no], S, Gm);
            v.push_back({Vertical(P.z[i]), S, pw[no], pa[no]});
        }
        return v;
    }
};

static void Passo2(const Opcoes &op) {
    using Clock = std::chrono::steady_clock;
    const auto t0 = Clock::now();
    const Zona &zo = Zonas()[op.zona - 1];
    std::printf("Kang et al. (2020), Passo 2: a coluna sem chuva; zona %d, modelo %s, ny = %d (celulas de %.3f m), "
                "dt = %g h\n", zo.id, op.tp ? "TP (bifasico, ar ativo)" : "SP (so a agua, ar passivo)", op.ny,
                H / op.ny, op.dt);
    Coluna col(zo, op);
    std::printf("  malha: %d elementos 2D, %lld equacoes\n", 4 * op.ny, (long long)col.m.CompMesh()->NEquations());
    RainTransient &chuva = *col.chuva;
    if (!chuva.InitialState()) {
        std::printf("  o estado inicial (geostatico) nao convergiu: %s\n", chuva.LastFailure().c_str());
        return;
    }
    TPZGeoMesh *g = col.m.GeoMesh();
    const TPZUnsaturatedHydraulics &h = col.Hidraulica();
    const REAL n = zo.L().thetaS;
    const Perfil &P = col.P;

    // a tabela de cada hora: pw e pa em 5 profundidades, a agua guardada e as conferencias
    const std::vector<REAL> zs = {0., 0.5, 1., 1.5, 2.};
    auto valor = [&](const std::vector<REAL> &v, REAL z) { // no perfil, o no mais perto de z
        size_t k = 0;
        for (size_t i = 1; i < P.z.size(); i++)
            if (std::fabs(P.z[i] - z) < std::fabs(P.z[k] - z)) k = i;
        return v[P.nos[k]];
    };
    const std::vector<REAL> pwI = chuva.NodalWaterPressure(), paI = chuva.NodalAirPressure();
    const REAL V0 = AguaGuardada(g, pwI, paI, h, n);
    std::printf("  estado inicial: p_w de %.2f a %.2f kPa, p_a de %.2f a %.2f kPa; agua guardada %.6f m3/m; "
                "|u| max %.2e m\n\n", chuva.PressureRange().first, chuva.PressureRange().second,
                chuva.AirPressureRange().first, chuva.AirPressureRange().second, V0, chuva.MaxDisplacement());
    std::printf("  %4s | %-39s | %-39s | %10s %9s %9s %6s %6s\n", "t(h)", "p_w (kPa) em z = 0 / 0,5 / 1 / 1,5 / 2 m",
                "p_a (kPa) em z = 0 / 0,5 / 1 / 1,5 / 2 m", "dV/V0", "simetria", "|u| max", "passos", "Newton");
    for (int hora = 0; hora <= int(duracao); hora++) {
        if (hora > 0) {
            chuva.SetIntensity(0.); // Passo 2: sem chuva
            if (!chuva.AdvanceTo(hora)) {
                std::printf("  o transiente parou em %.3f h: %s\n", chuva.Time(), chuva.LastFailure().c_str());
                return;
            }
        }
        const std::vector<REAL> pw = chuva.NodalWaterPressure(), pa = chuva.NodalAirPressure();
        // a simetria: a maior diferenca de p_w entre os nos esquerdo (2j) e direito (2j + 1) de uma mesma cota
        REAL sim = 0.;
        for (int j = 0; j <= op.ny; j++) sim = std::max(sim, std::fabs(pw[2 * j] - pw[2 * j + 1]));
        const REAL V = AguaGuardada(g, pw, pa, h, n);
        if (hora % 3 && hora != int(duracao)) continue; // imprime de 3 em 3 h e no fim
        std::printf("  %4d |", hora);
        for (REAL z : zs) std::printf(" %7.2f", valor(pw, z));
        std::printf(" |");
        for (REAL z : zs) std::printf(" %7.2f", valor(pa, z));
        std::printf(" | %10.2e %9.2e %9.2e %6lld %6lld\n", (V - V0) / V0, sim, chuva.MaxDisplacement(),
                    (long long)chuva.Steps(), (long long)chuva.NewtonIterations());
    }
    std::printf("\n  conferencias: dV/V0 perto de 0 (fronteiras fechadas; o resto e a compressibilidade da agua e a "
                "deformacao), simetria ~0 (1D), a agua desce para a base (p_w sobe embaixo e cai em cima)\n");
    std::printf("  tempo: %.1f s\n", std::chrono::duration<double>(Clock::now() - t0).count());
}

// ======================================================================================================== Passo 3
// A chuva hora a hora e o FS da Eq. 14 a cada hora
static void Passo3(const Opcoes &op) {
    using Clock = std::chrono::steady_clock;
    const auto t0 = Clock::now();
    const Zona &zo = Zonas()[op.zona - 1];
    const char *mod = op.tp ? "tp" : "sp";
    std::printf("Kang et al. (2020), Passo 3: a chuva da Fig. 2; zona %d, modelo %s, ny = %d, dt = %g h; angulos "
                "%.1f / %.1f / %.1f graus\n", zo.id, op.tp ? "TP (bifasico)" : "SP (so a agua)", op.ny, op.dt,
                zo.beta[0], zo.beta[1], zo.beta[2]);
    Coluna col(zo, op);
    if (op.normal)
        std::printf("  a coluna normal ao talude (Passo 7): beta = %.1f graus, espessura normal %.3f m, gravidade g cos(beta) "
                    "ao longo dela, chuva I cos(beta); %lld equacoes\n", col.beta * 180. / M_PI, H * std::cos(col.beta),
                    (long long)col.m.CompMesh()->NEquations());
    if (op.faixa)
        std::printf("  a faixa do FLAC (Passo 6): %g m na horizontal x %g m na vertical, beta = %.1f graus, %d x %d zonas, "
                    "%lld equacoes; o perfil e a vertical do meio (x = %g m)\n", Coluna::Lfaixa, H,
                    col.beta * 180. / M_PI, op.nx, op.ny, (long long)col.m.CompMesh()->NEquations(), Coluna::Lfaixa / 2);
    RainTransient &chuva = *col.chuva;
    if (!chuva.InitialState()) {
        std::printf("  o estado inicial (geostatico) nao convergiu: %s\n", chuva.LastFailure().c_str());
        return;
    }
    const REAL n = zo.L().thetaS, phi = zo.phi * M_PI / 180.;
    const std::string pre = op.saida + "_zona" + std::to_string(zo.id) + "_" + mod;
    std::ofstream fp(pre + "_perfis.csv"), ff(pre + "_fs.csv");
    fp << "t_h,z_m,pw_kPa,pa_kPa,Sw,succao_kPa\n";
    ff << "t_h,I_mmh,chuva_acum_mm,infiltrado_mm,face,FS_max,FS_medio,FS_min,zmin_medio_m,passos,newton,FS_beta\n";
    std::printf("  %4s %6s %8s %8s | %-39s | %-21s %6s %4s %6s\n", "t(h)", "I", "chuva", "infiltr",
                "p_w (kPa) em z = 0 / 0,5 / 1 / 1,5 / 2 m", "FS max / medio / min", "z_min", "face", "passos");
    REAL acum = 0.;
    for (int hora = 0; hora <= int(duracao); hora++) {
        if (hora > 0) {
            const REAL I = Chuva()[hora - 1]; // a chuva da hora anterior (mm/h)
            // na coluna normal, o topo e a superficie inclinada: a chuva vertical da I cos(beta) por metro dela
            chuva.SetIntensity(op.normal ? I * std::cos(col.beta) : I);
            acum += I;
            if (!chuva.AdvanceTo(hora)) {
                std::printf("  o transiente parou em %.3f h: %s\n", chuva.Time(), chuva.LastFailure().c_str());
                return;
            }
        }
        const std::vector<PontoPerfil> perf = col.PerfilAgora();
        REAL fs[3], zmin = 0.;
        for (int k = 0; k < 3; k++)
            fs[k] = FSInfinito(perf, zo.beta[k] * M_PI / 180., phi, zo.c, zo.gammaD, n, gammaW, nullptr,
                               k == 1 ? &zmin : nullptr);
        // o FS no angulo desta rodada: o da faixa (Passo 6) ou, na coluna, o angulo medio da zona
        const REAL fsBeta = (op.faixa || op.normal) ? FSInfinito(perf, col.beta, phi, zo.c, zo.gammaD, n, gammaW) : fs[1];
        const REAL infil = chuva.Infiltrated() / col.Largura() * 1000.; // m3/m -> mm
        for (const PontoPerfil &p : perf)
            fp << hora << "," << p.z << "," << p.pw << "," << p.pa << "," << p.Sw << "," << p.pa - p.pw << "\n";
        ff << hora << "," << (hora ? Chuva()[hora - 1] : 0.) << "," << acum << "," << infil << ","
           << chuva.FaceNodes() << "," << fs[0] << "," << fs[1] << "," << fs[2] << "," << zmin << ","
           << chuva.Steps() << "," << chuva.NewtonIterations() << "," << fsBeta << "\n";
        auto pz = [&](REAL z) { return perf[size_t(std::lround(z / H * op.ny))].pw; };
        std::printf("  %4d %6.1f %8.1f %8.1f | %7.2f %7.2f %7.2f %7.2f %7.2f | %6.3f %6.3f %6.3f %6.2f %4zu %6lld\n",
                    hora, hora ? Chuva()[hora - 1] : 0., acum, infil, pz(0.), pz(0.5), pz(1.), pz(1.5), pz(2.), fs[0],
                    fs[1], fs[2], zmin, chuva.FaceNodes(), (long long)chuva.Steps());
    }
    std::printf("  (face: nos do topo na face de percolacao, p_w = 0; 2 = a lamina d'agua: a chuva passa da "
                "capacidade de infiltracao e o excesso escorre)\n");
    std::printf("  CSV: %s_perfis.csv e %s_fs.csv; tempo: %.1f s\n", pre.c_str(), pre.c_str(),
                std::chrono::duration<double>(Clock::now() - t0).count());
}

// A malha geometrica em CSV: uma linha por elemento (o id do material ou da fronteira, o numero de vertices e as
// coordenadas x y de cada vertice), e o numero de equacoes da malha computacional no cabecalho
static void GravaMalha(Coluna &col, const std::string &arq) {
    std::ofstream f(arq);
    f << "# coluna do Kang2020: ny = " << col.op.ny << ", " << col.m.CompMesh()->NEquations()
      << " equacoes (u P2, p_w P1, p_a P1); ids: 1 solo, -1 base, -2 direita, -5 esquerda, -6 topo\n";
    f << "id,nv,x0,y0,x1,y1,x2,y2\n";
    for (TPZGeoEl *gel : col.m.GeoMesh()->ElementVec()) {
        if (!gel) continue;
        f << gel->MaterialId() << "," << gel->NCornerNodes();
        for (int i = 0; i < gel->NCornerNodes(); i++)
            f << "," << gel->NodePtr(i)->Coord(0) << "," << gel->NodePtr(i)->Coord(1);
        f << "\n";
    }
    std::printf("malha gravada em %s: %lld elementos, %lld equacoes\n", arq.c_str(),
                (long long)col.m.GeoMesh()->NElements(), (long long)col.m.CompMesh()->NEquations());
}

int main(int argc, char **argv) {
    const Opcoes op = LeOpcoes(argc, argv);
    if (!op.malha.empty()) {
        Coluna col(Zonas()[op.zona - 1], op);
        GravaMalha(col, op.malha);
        return 0;
    }
    if (op.passo == 1) Passo1();
    else if (op.passo == 2) Passo2(op);
    else Passo3(op);
    return 0;
}
