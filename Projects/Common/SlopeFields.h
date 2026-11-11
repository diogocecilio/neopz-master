// Campos aleatorios lognormais de c', phi' e k_v (Vargas Ceron et al. 2025, Secao 5.1.1): Karhunen-Loeve da covariancia
// exponencial separavel C = exp(-|dx|/Lx - |dy|/Ly) (Eq. 64) no retangulo que contem a malha (KLSeparavel.h, autopares em
// forma fechada), os mesmos modos para os tres campos, independentes entre si, e a transformacao lognormal das Eqs. 65-66.
// O campo e avaliado em qualquer posicao: as malhas de todos os niveis de refinamento veem a mesma realizacao
#ifndef SLOPEFIELDS_H
#define SLOPEFIELDS_H

#include "KLSeparavel.h"
#include "SlopeMesh.h"

#include <cstdint>

namespace slope {

// Estatistica dos campos; padrao: a configuracao de referencia do artigo (Tabela 2)
struct FieldStats {
    REAL c = 10., covC = 0.3;     // coesao (kPa)
    REAL phiDeg = 30., covPhi = 0.1; // atrito (graus)
    REAL covK = 0.6;              // k_v: o campo e o multiplicador k_v / media
    REAL Lx = 20., Ly = 2.;       // distancias de autocorrelacao (m)
    int M = 2000;                 // termos da expansao (o artigo usa 2000)
};

// Valores de uma realizacao num ponto
struct FieldValues {
    REAL c, phi, kmult; // kPa, rad, k_v / media
};

class RandomFields {
public:
    // Modos de KL no retangulo dos nos de m (calculados uma vez)
    RandomFields(const FieldStats &s, const FlatMesh &m);

    int M() const { return fKL.M(); }
    REAL VarianceFraction() const { return fKL.fracao; } // variancia retida, media no retangulo
    const FieldStats &Stats() const { return fStats; }

    // Sorteia a realizacao: tres sequencias normais padrao independentes (c', phi', k), com sementes derivadas de
    // (semente, realizacao, campo), como no SlopeRandomField (a realizacao n e a mesma la, com a mesma semente e retangulo)
    void Sample(uint64_t seed, int64_t realization);

    FieldValues At(REAL x, REAL y) const;

private:
    FieldStats fStats;
    klsep::KL2D fKL;
    std::vector<int> fUsedX, fUsedY;   // modos 1D usados por algum termo
    std::vector<int> fColX, fColY;     // modo 1D -> posicao em fUsedX / fUsedY
    std::vector<REAL> fXi[3];          // xi de cada campo (tamanho M)
    klsep::Lognormal fLc, fLphi, fLk;
};

} // namespace slope

#endif
