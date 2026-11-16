// Campos aleatorios lognormais pela expansao de Karhunen-Loeve analitica (KLSeparavel.h): covariancia exponencial
// separavel C = exp(-|dx|/Lx - |dy|/Ly) (Vargas Ceron et al. 2025, Eq. 64) num retangulo que contem o dominio, com os
// autopares em forma fechada, e a transformacao lognormal das Eqs. 65-66.
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
    // realizacao. O kernel e o mesmo (Lx, Ly do construtor); o sorteio, o deslocamento e o peso nao mudam
    void UseNystrom(const std::string &arquivo, int M);
    bool Nystrom() const { return !fNysLam.empty(); }
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
    // base de Nystrom (vazia: a analitica)
    std::vector<std::array<double, 3>> fNysNos;           // x, y, peso
    std::vector<double> fNysLam;                          // autovalores (fM)
    std::vector<std::vector<double>> fNysF;               // f_k nos nos: [no][k]
    double fNysFrac = 0.;                                  // sum lambda / area
    mutable std::vector<std::vector<double>> fNysA;       // a_j de cada campo na realizacao atual
    mutable bool fNysSujo = true;                          // os xi mudaram: refazer fNysA
};
