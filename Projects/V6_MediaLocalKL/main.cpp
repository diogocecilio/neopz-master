// V6_MediaLocalKL: executavel do caso completo (ver README.md). Uma execucao = um caso (um talude e o seu arco):
// constroi a KL de Comum/KLSeparavel.h no retangulo do talude (Fig. 9), calcula a reducao de variancia da media no arco
// pela Eq. 14 e pela KL truncada, a convergencia com o numero de termos e o Monte Carlo de N realizacoes.
//
// Variaveis de ambiente (padrao entre parenteses):
//   CASO       b40_h10 ... b40_h50, h50_b30 ... h50_b50, fig15 (b40_h10); "tabelas": so os arcos de Bishop e as
//              Eqs. 14/15 dos 10 casos, contra as Tabelas 5 a 8 (sem KL; segundos)
//   CONVENCAO  eq13 (L = delta / 2, a Eq. 13 do artigo) ou delta (L = delta) (eq13)
//   FRACAO     fracao de variancia da KL (0.90, a Eq. 29 do artigo);  MKL  limite de termos (400000)
//   N          realizacoes (10000);  SEMENTE (20210607);  NCONF  realizacoes conferidas por KL2D::Campo (20)
//   DELTA1, DELTA2  distancias de correlacao (3, 1 m);  VERBOSE (1)
// Saidas na pasta corrente: resumo.csv, convergencia.csv, sondas.csv, correlacao.csv (ou tabelas.csv) e, por ultimo,
// resultado.txt ("chave=valor ... situacao=ok seg=..."). Codigo de saida 0 em sucesso.
#include "problema.h"

#include <iostream>

int main() {
    try {
        const auto t0 = std::chrono::steady_clock::now();
        const v6::Config cf = v6::Config::DoAmbiente();
        std::ostream *log = cf.verbose ? &std::cout : nullptr;
        if (cf.caso == "tabelas") {
            const auto L = v6::CalculaTabelas(cf);
            const double seg = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
            v6::GravaTabelas(L, cf, seg);
            if (log) {
                std::cout << std::fixed << std::setprecision(4) << "V6 tabelas (convencao " << cf.convencao
                          << ", gamma calibrado = " << v6::GammaCalibrado() << " kN/m3)\n"
                          << "caso      arco    r(m)     T(m)   FS_B   FS_art  Gamma   s_cL  (art)   s_phiL (art)\n";
                for (const auto &l : L)
                    std::cout << std::setw(9) << std::left << l.caso.nome << std::right << std::setw(7)
                              << (l.caso.fig15 ? "fig15" : "bishop") << std::setw(8) << l.ac.arco.r << std::setw(8)
                              << l.ac.arco.T << std::setw(8) << l.ac.FSbishop << std::setw(8) << l.caso.FS
                              << std::setw(8) << l.G << std::setw(7) << l.G * v6::kSigmaC << std::setw(7)
                              << l.caso.scL << std::setw(8) << l.G * v6::kSigmaPhi << std::setw(7) << l.caso.sphiL
                              << "\n";
            }
            return 0;
        }
        const auto R = v6::RodaCaso(cf, log);
        v6::GravaCaso(R, cf);
        if (log) std::cout << "V6 " << cf.caso << ": fim (" << R.seg << " s)\n";
        return 0;
    } catch (const std::exception &e) {
        std::cerr << "V6_MediaLocalKL: erro: " << e.what() << "\n";
        return 1;
    }
}
