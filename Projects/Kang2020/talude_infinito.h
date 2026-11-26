// O fator de seguranca do talude infinito de Kang et al. (2020), Eqs. (14) e (15), sobre um perfil vertical:
//   FS(z) = { [W(z) cos^2(b) - (S_w p_w + S_g p_a)] tan(phi') + c' } / [ W(z) sin(b) cos(b) ]
//   W(z)  = int_0^z [ gamma_d + n S_w(z') gamma_w ] dz'
// com z a profundidade VERTICAL do plano de ruptura (paralelo a superficie), b o angulo do talude, S_g = 1 - S_w e o
// chi de Bishop igual a S_w (Eq. 12). O FS do talude e o minimo de FS(z) ao longo do perfil.
#pragma once

#include "pzreal.h"

#include <cmath>
#include <limits>
#include <vector>

namespace kang {

/// Um ponto do perfil: profundidade vertical z (m, 0 na superficie), saturacao e pressoes da agua e do ar (kPa)
struct PontoPerfil {
    REAL z, Sw, pw, pa;
};

/// O FS minimo do perfil (pontos em z crescente, o primeiro na superficie). beta e phi em radianos. Opcional: o FS
/// em cada ponto (fsz; infinito na superficie, onde W = 0) e a profundidade do minimo (zmin)
inline REAL FSInfinito(const std::vector<PontoPerfil> &perfil, REAL beta, REAL phi, REAL c, REAL gammaD, REAL n,
                       REAL gammaW, std::vector<REAL> *fsz = nullptr, REAL *zmin = nullptr) {
    const REAL cb = std::cos(beta), sb = std::sin(beta), tf = std::tan(phi);
    REAL W = 0., fsMin = std::numeric_limits<REAL>::infinity();
    if (fsz) fsz->assign(perfil.size(), std::numeric_limits<REAL>::infinity());
    for (size_t i = 1; i < perfil.size(); i++) {
        const PontoPerfil &a = perfil[i - 1], &b = perfil[i];
        // Eq. (15) pela regra do trapezio: o peso da coluna acima de z, por unidade de area horizontal
        W += 0.5 * ((gammaD + n * a.Sw * gammaW) + (gammaD + n * b.Sw * gammaW)) * (b.z - a.z);
        const REAL pporo = b.Sw * b.pw + (1. - b.Sw) * b.pa; // a pressao de poros de Bishop com chi = S_w (Eq. 12)
        const REAL fs = ((W * cb * cb - pporo) * tf + c) / (W * sb * cb);
        if (fsz) (*fsz)[i] = fs;
        if (fs < fsMin) {
            fsMin = fs;
            if (zmin) *zmin = b.z;
        }
    }
    return fsMin;
}

} // namespace kang
