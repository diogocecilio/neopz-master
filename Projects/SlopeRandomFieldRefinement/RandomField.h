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

    int M() const { return fKL.M(); }
    int NCampos() const { return fLog.size(); }
    double VarianceFraction() const { return fKL.fracao; } // variancia retida, media no retangulo

    // Sorteia a realizacao: uma sequencia normal padrao por campo, com semente derivada de (semente, realizacao, campo).
    // Com 3 campos, a mesma formula do SlopeFields do Common: a realizacao n e a mesma de la
    void Sample(uint64_t semente, int64_t realizacao);

    // Valores de todos os campos no ponto (na ordem dos campos do construtor)
    std::vector<double> At(double x, double y) const;

    // ---- amostragem por importancia: deslocamento da media no espaco dos xi
    // Direcao unitaria d de um conjunto de pontos {x, y, peso}: d ~ a, com a_k = sqrt(lambda_k) f_k(x) g_k(y) na media
    // ponderada dos pontos, o gradiente da media do campo gaussiano Z nesses pontos em relacao aos xi (a mesma nos
    // campos, que usam os mesmos modos). Vazia se |a| = 0
    std::vector<double> Direction(const std::vector<std::array<double, 3>> &pontos) const;
    // Desloca a media de cada campo c: xi_c ~ N(s_c d, I) no Sample. s vazio volta ao sorteio sem deslocamento
    void SetShift(const std::vector<double> &d, const std::vector<double> &s);
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
};
