// Transiente do rebaixamento do reservatorio com esqueleto rigido (u = 0) e ar passivo (p_a = 0): so a agua percola,
// por Euler implicito, com o TPZUnsatAnalysis do neopz (como o RainTransient). O nivel do reservatorio desce em linha
// reta de nivelInicial ate nivelFinal em tRebaixamento horas e depois fica parado. Contornos:
//  - lateral do macico (idMacico): p_w hidrostatica com o nivel nivelMacico (fixo) abaixo dele; acima, sem fluxo;
//  - lado do reservatorio (idsReservatorio): abaixo do nivel do reservatorio, p_w = gamma_w (nivel - y); acima, face de
//    percolacao: p_w = 0 onde a agua sai do solo, sem fluxo onde nao sai;
//  - base e o resto: sem fluxo.
// Estado inicial: hidrostatico com o nivel do reservatorio (o macico em equilibrio com o reservatorio cheio)
//
// Contorno do Ceron (contornoCeron, construtor com DrawdownCeronIds): os contornos das Eqs. 19-21 de Vargas Ceron et
// al. (2025), em excesso u = p_w - gamma_w (nivelInicial - y) (nivelInicial = a crista), todos Dirichlet, sem face de
// percolacao:
//  - crista: u = 0 (p_w = 0); lateral esquerda e base: u = 0 (p_w = gamma_w (nivelInicial - y));
//  - face e terreno do pe: p_w = gamma_w max(nivel - y, 0) com o nivel do reservatorio do instante (no pe, com
//    h_w <= 1, e gamma_w (nivel - y), isto e u = -gamma_w h_w);
//  - lateral direita (e o resto): sem fluxo.
#pragma once

#include "TPZUnsatAnalysis.h"

#include <map>
#include <set>
#include <vector>

struct DrawdownCeronIds { // ids de contorno para o contorno do Ceron
    int crista = 0, face = 0, pe = 0, base = 0, esquerda = 0;
};

struct DrawdownData {
    REAL gammaw = 9.81;
    REAL nivelInicial = 0., nivelFinal = 0.; // cotas do reservatorio no inicio e no fim do rebaixamento (m)
    REAL nivelMacico = 0.;                    // cota do NA na lateral do macico (m)
    REAL tRebaixamento = 24.;                 // duracao do rebaixamento (h)
    REAL dt = 0.1;                            // passo de tempo (h)
    bool contornoCeron = false;               // contornos do Ceron (posto pelo construtor com DrawdownCeronIds)
    REAL Nivel(REAL tHoras) const             // cota do reservatorio no instante t (h)
    {
        const REAL f = tRebaixamento > 0. ? std::min(tHoras / tRebaixamento, REAL(1)) : 1.;
        return nivelInicial + f * (nivelFinal - nivelInicial);
    }
};

class DrawdownTransient {
public:
    DrawdownTransient(TPZUnsatAnalysis &an, const DrawdownData &d, int idMacico, const std::set<int> &idsReservatorio);
    // contorno do Ceron (liga contornoCeron)
    DrawdownTransient(TPZUnsatAnalysis &an, const DrawdownData &d, const DrawdownCeronIds &ids);

    void InitialState();          // t = 0: p_w hidrostatica com o nivel inicial do reservatorio
    bool AdvanceTo(REAL tHoras);  // passos de dt ate t; false se um passo nao convergiu
    REAL Time() const { return fTime / 3600.; }
    REAL Level() const { return fData.Nivel(Time()); }
    std::vector<REAL> NodalWaterPressure() const; // p_w por no da malha geometrica

private:
    bool Step(REAL dt);      // um passo, corrigindo a face de percolacao ate ela se repetir
    bool StepCeron(REAL dt); // um passo com os contornos do Ceron (todos Dirichlet)
    void FixedBase();        // u = 0 e p_a = 0

    TPZUnsatAnalysis &fAn;
    DrawdownData fData;
    std::map<int64_t, REAL> fFixed;    // sempre impostas: u = 0, p_a = 0 e a lateral do macico
    std::set<int64_t> fReservatorio;  // nos do lado do reservatorio
    std::set<int64_t> fFace;           // nos na face de percolacao no ultimo passo aceito
    std::set<int64_t> fFaceCeron;      // contorno do Ceron: nos da face e do terreno do pe (p_w do nivel do instante)
    REAL fTime = 0.;                   // s
};
