// V3_SuperficieLivre: o caso completo de um exemplo de Sharma et al. (2021) ate o regime permanente (ver README.md).
// Configuracao por variaveis de ambiente (padrao entre parenteses):
//   CASO      retangular (4.1), inclinada (4.2) ou dois_blocos (4.3)            (retangular)
//   NX, NY    divisoes da malha estruturada                                     (retangular 32 x 48; inclinada 28 x 40;
//                                                                                 dois_blocos 40 x 80)
//   RX        razao ultimo/primeiro intervalo em x (< 1 refina junto a jusante)  (1)
//   RET       curva de retencao: vg (van Genuchten) ou gardner                  (vg)
//   AG, BVG   vg: S_e = [1 + (AG p_c)^BVG]^(1/BVG - 1); gardner: S_e = exp(-AG p_c) (AG em 1/kPa)   (1 e 3)
//   POT       k_rw = S_e^POT                                                    (3)
//   KRW_MIN   limite inferior de k_rw                                           (1e-9)
//   POROS     porosidade (so o pseudo-transiente)                               (0.3)
//   DT0, DT_FATOR, DT_MAX, T_MAX  passos do pseudo-transiente (s)               (1, 1.5, 1e5, 1e7)
//   TOL_EST   regime permanente: |q_entra - q_sai| / q_entra                    (1e-4)
//   CI        estado inicial: hidrostatica ou linear                            (hidrostatica)
//   TOPO_FACE 1: a crista tambem pode ser face de percolacao                    (1)
//   VERBOSE   0, 1 (cada passo) ou 2 (cada iteracao de Newton)                  (1)
// Saidas na pasta corrente: passos.csv (cada passo aceito: vazoes e face de percolacao), face.csv (nos da face a
// jusante no estado final), linha_freatica.csv (a nossa linha p_w = 0 e a digitalizada do artigo), pw.csv (p_w nos
// vertices), triangulos.csv, resumo.csv e resultado.txt (gravado por ultimo; um .tmp renomeado).
#include "problema.h"

#include <cstdio>
#include <cstdlib>

static double envd(const char *k, double d) { const char *v = std::getenv(k); return v && *v ? std::atof(v) : d; }
static int envi(const char *k, int d) { const char *v = std::getenv(k); return v && *v ? std::atoi(v) : d; }
static std::string envs(const char *k, const std::string &d) { const char *v = std::getenv(k); return v && *v ? v : d; }

int main() {
    const auto inicio = v3::Relogio::now();
    const v3::Caso c = v3::CasoPorNome(envs("CASO", "retangular"));
    if (c.nome.empty()) return 2;
    v3::Numerica num;
    // malhas padrao do caso completo: espacamento de 0,5 m (retangular), 0,125 m na vertical (inclinada e dois blocos)
    if (c.nome == "inclinada") num.nx = 28, num.ny = 40;
    if (c.nome == "dois_blocos") num.nx = 40, num.ny = 80;
    num.nx = envi("NX", num.nx), num.ny = envi("NY", num.ny);
    num.rx = envd("RX", num.rx);
    num.ret = envs("RET", num.ret), num.bVG = envd("BVG", num.bVG);
    num.aG = envd("AG", num.aG), num.pot = envd("POT", num.pot), num.krwMin = envd("KRW_MIN", num.krwMin);
    num.n = envd("POROS", num.n);
    num.dt0 = envd("DT0", num.dt0), num.fator = envd("DT_FATOR", num.fator), num.dtMax = envd("DT_MAX", num.dtMax);
    num.tMax = envd("T_MAX", num.tMax), num.tolEst = envd("TOL_EST", num.tolEst);
    num.ci = envs("CI", num.ci);
    num.topoFace = envi("TOPO_FACE", 1) != 0;
    num.verbose = envi("VERBOSE", 1);
    if (num.nx < 2 || num.ny < 2 || !(num.aG > 0.) || !(num.dt0 > 0.) || !(num.fator >= 1.)) {
        std::cerr << "parametros invalidos (NX, NY, AG, DT0, DT_FATOR)\n";
        return 2;
    }
    std::cout << "V3_SuperficieLivre: caso " << c.nome << " (Sharma et al. 2021, Sec. " << c.secao << "), malha " << num.nx
              << " x " << num.ny << ", retencao " << num.ret << " a = " << num.aG << " 1/kPa"
              << (num.ret == "vg" ? " b = " + std::to_string(num.bVG) : std::string()) << ", k_rw = S_e^" << num.pot
              << std::endl;

    std::ofstream fpas("passos.csv");
    fpas << "passo,t_s,dt_s,q_esq_m3s,q_jus_m3s,q_face_m3s,q_topo_m3s,dVdt_m3s,desequilibrio,face_nos,y_saida_m,"
            "tentativas_face,newton,seg\n"
         << std::setprecision(10);
    const v3::Resultado r = v3::Executa(c, num, &fpas);
    fpas.close();

    // face a jusante no estado final
    {
        std::ofstream f("face.csv");
        f << "x_m,y_m,pw_kPa,q_saida_m_s,tipo\n" << std::setprecision(10);
        for (const auto &n : r.face) f << n.x << ',' << n.y << ',' << n.pw << ',' << n.q << ',' << n.tipo << '\n';
    }
    // linha freatica: a nossa (p_w = 0 na vertical x) e a digitalizada do artigo, quando houver
    std::vector<std::array<REAL, 2>> ref;
    if (!c.csvLinha.empty()) ref = v3::LeLinha(std::string(VERIF_DIR) + "/dados/" + c.csvLinha);
    {
        std::ofstream f("linha_freatica.csv");
        f << "x_m,y_nosso_m,y_artigo_m\n" << std::setprecision(8);
        const REAL xmax = std::max(c.xB, c.xD);
        const int nAm = 160;
        for (int i = 0; i <= nAm; ++i) {
            const REAL x = xmax * i / nAm;
            const REAL y = r.Freatica(std::min(x, xmax * (1. - 1.e-9)));
            f << x << ',' << y << ',' << (ref.empty() ? -1. : v3::Interpola(ref, x)) << '\n';
        }
    }
    {
        std::ofstream f("pw.csv");
        f << "no,x_m,y_m,pw_kPa\n" << std::setprecision(10);
        for (size_t i = 0; i < r.nosPw.size(); ++i)
            f << i << ',' << r.nosPw[i][0] << ',' << r.nosPw[i][1] << ',' << r.nosPw[i][2] << '\n';
        std::ofstream t("triangulos.csv");
        t << "n0,n1,n2\n";
        for (const auto &tr : r.tri) t << tr[0] << ',' << tr[1] << ',' << tr[2] << '\n';
    }
    // desvio medio da linha freatica em relacao a digitalizada (nos pontos digitalizados)
    REAL desvioMax = -1., desvioMed = 0.;
    int nDesv = 0;
    for (const auto &p : ref) {
        if (p[0] <= 0. || p[0] >= std::max(c.xB, c.xD) - 1.e-9) continue;
        const REAL y = r.Freatica(p[0]);
        if (y < 0.) continue;
        desvioMax = std::max(desvioMax, std::fabs(y - p[1]));
        desvioMed += std::fabs(y - p[1]);
        nDesv++;
    }
    if (nDesv) desvioMed /= nDesv;
    const REAL qSai = r.ultimo.qJus + r.ultimo.qFace + r.ultimo.qTopo;
    const double seg = v3::Segundos(inicio);
    const std::string situacao = r.ok ? "ok" : "falhou";

    std::ostringstream lin;
    lin << std::setprecision(8) << "caso=" << c.nome << " nx=" << num.nx << " ny=" << num.ny << " ret=" << num.ret
        << " aG=" << num.aG << " bVG=" << num.bVG
        << " pot=" << num.pot << " nos=" << r.nos << " passos=" << r.passos << " cortes=" << r.cortes
        << " newton=" << r.newton << " t_s=" << r.t << " y_saida=" << r.ySaida << " y_acima=" << r.yAcima
        << " y_saida_est=" << r.ySaidaEst() << " y_saida_artigo=" << c.ySaidaArtigo << " q_esq=" << r.ultimo.qEsq
        << " q_sai=" << qSai << " q_face=" << r.ultimo.qFace << " q_topo=" << r.ultimo.qTopo << " q_DC=" << r.qDC
        << " q_DC_corr=" << r.qDCcorr << " q_franja_max=" << r.qFranjaMax << " desequilibrio=" << r.ultimo.est
        << " y_interface=" << r.yInterface << " y_interface_pw0=" << r.yInterface0
        << " desvio_linha_med=" << desvioMed << " desvio_linha_max=" << desvioMax << " situacao=" << situacao
        << " seg=" << seg;
    {
        std::ofstream f("resumo.csv");
        f << "caso,nx,ny,aG,pot,nos,passos,newton,t_s,y_saida_m,y_acima_m,y_saida_artigo_m,q_esq_m3s,q_sai_m3s,q_DC_m3s,"
             "q_DC_corr_m3s,q_franja_max_m3s,desequilibrio,y_interface_m,y_interface_pw0_m,desvio_linha_med_m,"
             "desvio_linha_max_m,situacao,seg\n"
          << std::setprecision(10) << c.nome << ',' << num.nx << ',' << num.ny << ',' << num.aG << ',' << num.pot << ','
          << r.nos << ',' << r.passos << ',' << r.newton << ',' << r.t << ',' << r.ySaida << ',' << r.yAcima << ','
          << c.ySaidaArtigo << ',' << r.ultimo.qEsq << ',' << qSai << ',' << r.qDC << ',' << r.qDCcorr << ','
          << r.qFranjaMax << ',' << r.ultimo.est << ',' << r.yInterface << ',' << r.yInterface0 << ',' << desvioMed
          << ',' << desvioMax << ',' << situacao << ','
          << seg << '\n';
    }
    {
        std::ofstream f("resultado.txt.tmp");
        f << lin.str() << '\n';
        f.close();
        std::rename("resultado.txt.tmp", "resultado.txt");
    }
    std::cout << lin.str() << std::endl;
    if (!r.ok) std::cerr << "V3: " << r.msg << std::endl;
    return r.ok ? 0 : 1;
}
