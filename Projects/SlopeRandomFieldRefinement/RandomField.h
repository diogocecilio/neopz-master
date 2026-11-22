// Campos aleatorios lognormais pela expansao de Karhunen-Loeve analitica (KLSeparavel.h): covariancia exponencial
// separavel C = exp(-|dx|/Lx - |dy|/Ly) (Vargas Ceron et al. 2025, Eq. 64) num retangulo que contem o dominio, com os
// autopares em forma fechada, e a transformacao lognormal das Eqs. 65-66. Ou pela base de Nystrom do proprio dominio
// (UseNystrom: Projects/Cho2010/nystrom_base.py), com o nucleo exponencial separavel ou o gaussiano (girado ou nao).
// Varios campos independentes (ex.: c', phi' e k) usam os mesmos modos e cada um tem as suas variaveis normais xi.
// O campo e avaliado em qualquer ponto (x, y): malhas diferentes (a do transiente e a do fator) veem a mesma realizacao.
// Generico: nada de talude aqui; o problema diz o retangulo, as distancias de autocorrelacao e a media e o CoV de cada
// campo
#pragma once

#include "KLSeparavel.h"

#include <array>
#include <cstdint>
#include <string>
#include <vector>

class RandomField {
public:
    // media e coeficiente de variacao de um campo lognormal
    struct Stats {
        double media, cov;
    };

    // Modos de KL do retangulo [x0, x0 + Dx] x [y0, y0 + Dy] com as distancias de autocorrelacao Lx, Ly e M termos
    // (calculados uma vez); campos: media e CoV de cada campo
    RandomField(double x0, double y0, double Dx, double Dy, double Lx, double Ly, int M, const std::vector<Stats> &campos);

    int M() const { return fM; }
    int NCampos() const { return fLog.size(); }
    double VarianceFraction() const { return Nystrom() ? fNysFrac : fKL.fracao; } // variancia retida (media no dominio)

    // Base de Nystrom no lugar da analitica (arquivo do Projects/Cho2010/nystrom_base.py: os nos de quadratura x, y, w, os
    // autovalores e f_k nos nos): os M primeiros modos (M <= os do arquivo). O campo num ponto qualquer sai da interpolacao
    // de Nystrom, Z(x) = sum_j C(x, x_j) a_j, com a_j = w_j sum_k xi_k f_k(x_j) / sqrt(lambda_k) calculado uma vez por
    // realizacao; o sorteio, o deslocamento e o peso nao mudam. O nucleo C vem do cabecalho do arquivo (nucleo=exp ou
    // gauss, lh, lv, angulo); sem ele (bases antigas), o exponencial com os Lx, Ly do construtor.
    // rapido = true: com o nucleo exponencial e os nos numa grade regular (a do nystrom_base.py), a soma sobre os nos vira
    // quatro recorrencias na grade (exp(-|dx|/L) = produto de fatores exp(-dx/L) de celula em celula), calculadas uma vez
    // por realizacao: o campo num ponto custa O(1) em vez de O(n). E a mesma soma, em outra ordem (difere so no
    // arredondamento, ~1e-12); false (o padrao) deixa as contas de antes, bit a bit
    void UseNystrom(const std::string &arquivo, int M, bool rapido = false);
    bool Nystrom() const { return !fNysLam.empty(); }
    bool FastNystrom() const { return fGrade; }
    // o nucleo da base de Nystrom, por extenso (ex.: "exp lh=20 lv=2")
    std::string KernelDescription() const;
    double KernelLx() const { return fLx; }
    double KernelLy() const { return fLy; }
    bool GaussianKernel() const { return fNucleo == EGauss; }
    // soma_p a_p sqrt(lambda_k) f_k(x_p) (= soma_p a_p dZ/dxi_k no ponto p) para k < K, para cada vetor de pesos a
    // (pesos[i][p]): a projecao nos modos de uma sensibilidade dada nos pontos (o gradiente modal). Na base de Nystrom,
    // (1/sqrt(lambda_k)) sum_j w_j f_k(x_j) b_j com b_j = sum_p a_p C(x_p, x_j): O(P n + n K) em vez dos O(P n K) do
    // ModeValues ponto a ponto
    std::vector<std::vector<double>> ProjectModes(const std::vector<std::array<double, 2>> &pontos,
                                                  const std::vector<std::vector<double>> &pesos, int K) const;
    // sqrt(lambda_k) f_k(x, y) dos K primeiros termos (K <= M; a mesma base em todos os campos): dZ/dxi_k no ponto
    std::vector<double> ModeValues(double x, double y, int K) const;
    double SigmaLn(int c) const { return fLog[c].sigma_ln; }
    int CorrA() const { return fCorrA; }
    int CorrB() const { return fCorrB; }
    double Rho() const { return fRho; }

    // Sorteia a realizacao: uma sequencia normal padrao por campo, com semente derivada de (semente, realizacao, campo).
    // Com 3 campos, a mesma formula do SlopeFields do Common: a realizacao n e a mesma de la
    void Sample(uint64_t semente, int64_t realizacao);

    // Valores de todos os campos no ponto (na ordem dos campos do construtor)
    std::vector<double> At(double x, double y) const;

    // Correlacao cruzada rho entre os campos gaussianos a e b (antes da transformacao lognormal; como em Cho 2010, Eq. 11):
    // o campo b passa a ser rho Z_a + sqrt(1 - rho^2) Z_b, com a mesma covariancia espacial (os campos usam os mesmos
    // modos). Os xi continuam independentes: o sorteio, o deslocamento e o peso da amostragem por importancia nao mudam
    void SetCrossCorrelation(int a, int b, double rho);

    // ---- amostragem por importancia: deslocamento da media no espaco dos xi
    // Direcao unitaria d de um conjunto de pontos {x, y, peso}: d ~ a, com a_k = sqrt(lambda_k) f_k(x) g_k(y) na media
    // ponderada dos pontos, o gradiente da media do campo gaussiano Z nesses pontos em relacao aos xi (a mesma nos
    // campos, que usam os mesmos modos). Vazia se |a| = 0
    std::vector<double> Direction(const std::vector<std::array<double, 3>> &pontos) const;
    // Desloca a media de cada campo c: xi_c ~ N(s_c d, I) no Sample. s vazio volta ao sorteio sem deslocamento
    void SetShift(const std::vector<double> &d, const std::vector<double> &s);
    // O mesmo com o deslocamento de cada termo de cada campo dado (mu[c][k]; termos que faltam ficam em 0)
    void SetShift(const std::vector<std::vector<double>> &mu);
    // Os xi da ultima realizacao (com o deslocamento), por campo
    const std::vector<std::vector<double>> &Xi() const { return fXi; }
    // Os xi de todos os campos dados, sem sorteio (ex.: um ponto do FORM): ate M valores por campo, o resto 0
    void SetXi(const std::vector<std::vector<double>> &xi);
    // ln do peso da ultima realizacao, ln(phi(xi) / phi(xi - mu)) = sum_c (-xi_c . mu_c + |mu_c|^2 / 2): com ele o
    // estimador da probabilidade e sem vies (0 sem deslocamento)
    double LogWeight() const;
    // Projecao xi_c . d de cada campo na ultima realizacao
    std::vector<double> Projection(const std::vector<double> &d) const;

private:
    klsep::KL2D fKL;
    std::vector<int> fUsadosX, fUsadosY; // modos 1D usados por algum termo
    std::vector<int> fColX, fColY;       // modo 1D -> posicao em fUsadosX / fUsadosY
    std::vector<klsep::Lognormal> fLog;  // transformacao de cada campo
    std::vector<std::vector<double>> fXi; // xi de cada campo (tamanho M)
    std::vector<std::vector<double>> fMu; // deslocamento da media de cada campo (vazio: sem deslocamento)
    int fCorrA = -1, fCorrB = -1;         // correlacao cruzada entre os campos gaussianos fCorrA e fCorrB
    double fRho = 0.;
    int fM = 0;                           // termos usados
    double fLx = 1., fLy = 1.;            // distancias de autocorrelacao (o kernel do Nystrom)
    enum ENucleo { EExp = 0, EGauss = 1 };
    ENucleo fNucleo = EExp;               // nucleo da base de Nystrom
    double fAngulo = 0.;                  // gauss: rotacao dos eixos (rad)
    double Kernel(double dx, double dy) const; // C(dx, dy) do nucleo da base
    // grade regular dos nos (rapido): fGNx x fGNy celulas, centro da primeira (fGx0, fGy0), passo fGdx x fGdy; celula de
    // cada no; e, por campo, as quatro somas acumuladas da realizacao (frente/tras em y e depois em x)
    bool fGrade = false;
    int fGNx = 0, fGNy = 0;
    double fGx0 = 0., fGy0 = 0., fGdx = 1., fGdy = 1.;
    std::vector<int64_t> fGCelula;
    mutable std::vector<std::array<std::vector<double>, 4>> fGQ;
    void PrepareGrid() const;              // as somas da realizacao (fNysA ja feito)
    double GridValue(int c, double x, double y) const; // Z do campo c no ponto pelas somas
    // base de Nystrom (vazia: a analitica)
    std::vector<std::array<double, 3>> fNysNos;           // x, y, peso
    std::vector<double> fNysLam;                          // autovalores (fM)
    std::vector<std::vector<double>> fNysF;               // f_k nos nos: [no][k]
    double fNysFrac = 0.;                                  // sum lambda / area
    mutable std::vector<std::vector<double>> fNysA;       // a_j de cada campo na realizacao atual
    mutable bool fNysSujo = true;                          // os xi mudaram: refazer fNysA
};
