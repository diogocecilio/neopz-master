// Caminho de equilibrio do talude com o filtro de equacoes: Newton, arc-length (de Souza Neto et al. 2008, Box 4.4) e
// continuacao com bisseccao, e os fatores GIM (aumento do peso) e SRM (reducao da resistencia) por cada um deles
#ifndef SLOPESOLVER_H
#define SLOPESOLVER_H

#include "TPZLinearAnalysis.h"

#include <chrono>
#include <cstdio>
#include <functional>
#include <iostream>
#include <string>
#include <vector>

namespace slope {

using Clock = std::chrono::steady_clock;
inline double Seconds(Clock::time_point t0) { return std::chrono::duration<double>(Clock::now() - t0).count(); }

// Silencia o std::cout no escopo (a montagem serial do NeoPZ imprime '*' a cada 1000 elementos)
struct QuietCout {
    std::streambuf *old;
    QuietCout() : old(std::cout.rdbuf(nullptr)) {}
    ~QuietCout() { std::cout.rdbuf(old); }
};

// Linha formatada como printf
template <class... A>
void Print(const char *fmt, A... a)
{
    char buf[256];
    std::snprintf(buf, sizeof buf, fmt, a...);
    std::cout << buf << std::endl;
}

// O que o caminho usa de cada material: a analise, as equacoes livres (as outras ficam em x0, ou em lambda x0 com
// scaleFixed) e como mudar o peso (lambda) e a resistencia (c/F e tan(phi)/F)
struct PathProblem {
    TPZLinearAnalysis *an = nullptr;
    std::vector<int64_t> free, fixed;
    TPZFMatrix<STATE> x0;     // valores das fixas com lambda = 1 (pressoes congeladas; 0 no seco)
    bool scaleFixed = false;  // fixas = lambda x0 (u-p); senao o material escala as pressoes
    REAL lambda = 1., reduction = 1.;
    std::function<void(REAL)> setLambda, setReduction;
    std::function<void()> accept; // memoria dos pontos <- solucao carregada na malha
    void SetLambda(REAL l) { lambda = l, setLambda(l); }
    void SetReduction(REAL F) { reduction = F, setReduction(F); }
    REAL Dot(const TPZFMatrix<STATE> &a, const TPZFMatrix<STATE> &b) const;
};

void SetFreeEquations(TPZLinearAnalysis &an, const std::vector<int64_t> &free); // filtro da matriz
TPZFMatrix<STATE> SolveK(TPZLinearAnalysis &an, const TPZFMatrix<STATE> &b);   // K y = b, fatoracao corrente
bool SolveFree(TPZLinearAnalysis &an, const std::vector<int64_t> &free, TPZFMatrix<STATE> &x); // problema linear

enum class FactorMethod { GIMArcLength, GIMBisection, SRMArcLength, SRMBisection };
const char *Name(FactorMethod m);

struct FactorResult {
    REAL value = 0.; // Gamma (GIM) ou FS (SRM)
    int steps = 0;
    std::string end; // por que o caminho terminou
};

// Fator de seguranca a partir do estado sem tensoes (x0 nas fixas); imprime uma linha por passo
FactorResult ComputeFactor(PathProblem &p, int64_t neq, FactorMethod method);

} // namespace slope

#endif
