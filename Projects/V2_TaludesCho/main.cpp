// V2_TaludesCho: executavel do contrato de Monte Carlo (uma realizacao por execucao, na pasta corrente; ver README.md).
// Variaveis de ambiente principais: CASO=cho1|cho2, REAL=<n>, CAMPOS=0|1 (0: caso deterministico com as medias),
// SEMENTE, NUCLEO=gauss|exp, COV_C, COV_PHI, LX, LY, RHO, KL_FRAC, KL_M, H_MALHA, FS_TOL, MAXIT_FROZEN, CAMPOS_CSV,
// FS_METODO=bissecao|arco (o SRM por bissecao ou por arc-length).
// Saidas: fs.csv, fs_tentativas.csv, resumo.csv e, por ultimo, resultado.txt. Codigo de saida 0 em sucesso.
#include "problema.h"

#include <iostream>

int main() {
    v2::Verboso() = v2::EnvI("VERBOSO", 1);
    try {
        return v2::RodaExecucao(v2::Config::DoAmbiente(), ".");
    } catch (const std::exception &e) {
        std::cerr << "V2_TaludesCho: " << e.what() << std::endl;
        return 2;
    }
}
