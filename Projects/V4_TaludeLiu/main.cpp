// V4_TaludeLiu: executavel do contrato de Monte Carlo (uma realizacao por execucao, na pasta corrente; ver README.md):
// o talude de Liu et al. (2024) sob chuva permanente q/k_s = 0.01, regime permanente por um transiente longo (variante
// iii) e FS por reducao de resistencia. Variaveis de ambiente principais: CASO=det|s1|s2|s3 (det: as medias; s1: campo
// de k_s; s2: campos de c' e phi'; s3: os tres), ANGULO=0|90 (acamamento horizontal ou vertical: troca Lx e Ly),
// REAL=<n>, CAMPOS=0|1 (0: as medias), SEMENTE, COV_C, COV_PHI, COV_K, RHO, LX, LY, KL_M, MALHA, P_KR, Q_KS,
// PW_PERMANENTE (le o permanente de um pw.csv quando k_s e constante), VERBOSO. Saidas: passos.csv, pw.csv, secoes.csv,
// fs_tentativas.csv, fs.csv, resumo.csv e, por ultimo, resultado.txt. Codigo de saida 0 em sucesso.
#include "problema.h"

#include <iostream>

int main() {
    v4::Verboso() = v4::EnvI("VERBOSO", 1);
    try {
        const v4::Config cf = v4::Config::DoAmbiente();
        TPZUnsatAnalysis::SetUseSkyline(v4::EnvS("SOLVER", "pardiso") == "skyline");
        return v4::RodaExecucao(cf, ".").codigo;
    } catch (const std::exception &e) {
        std::cerr << "V4_TaludeLiu: " << e.what() << std::endl;
        return 2;
    }
}
