// Arquivos da animacao (o formato que o animacao.py le): a malha, os campos do transiente em cada instante gravado e,
// em cada instante com fator de seguranca, o estado do fim do GIM (deslocamento nos vertices; sqrt(J2(eps_p)),
// hardening e tipo de resposta nos pontos de integracao). Malha: os triangulos do dominio do modelo do fator (a malha do
// transiente tem que ser uma copia dela, com a mesma numeracao de nos). Arquivos, na pasta dada:
//   vertices.csv (no,x,y), triangles.csv (n0,n1,n2), pontos.csv (tri,x,y,peso), snaps_pw/pa/sw.csv (t e um valor por
//   vertice), fs.csv (t_h,FS), srm_nos.csv (t_h,no,ux,uy,ux_geo,uy_geo), srm_pontos.csv (t_h,ponto,raizJ2,hardening,tipo),
//   nivel.csv (t_h,nivel do reservatorio) e summary.txt
// Rebaixamento (main.cpp, animacao_rebaixamento.py): Snapshot com o nivel, Summary e FactorOfSafety(t, fs) do proprio
// modelo. Chuva (rahardjo.cpp, LearningMohrCoulomb2/animacao.py): o modelo e o do transiente, Snapshot sem nivel,
// RainSummary e FactorOfSafety(t, fs, fator), que leva o estado do fator, calculado em outra malha, para esta
#pragma once

#include "PointLocator.h"
#include "UnsatModel.h"

#include "Plasticity/TPZUnsaturatedHydraulics.h"
#include "TPZTensor.h"

#include "pzgeoel.h"
#include "pzintel.h"
#include "pzmultiphysicselement.h"
#include "pzquad.h"

#include <array>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <map>
#include <memory>
#include <string>
#include <vector>

template <class TMat>
class AnimationFiles {
public:
    // pasta: criada se nao existir; model: o modelo do fator (malha computacional ja criada)
    AnimationFiles(const std::string &pasta, UnsatModel<TMat> &model) : fPasta(pasta), fModel(model)
    {
        std::filesystem::create_directories(pasta);
        TPZGeoMesh *g = model.GeoMesh();
        std::ofstream v(Arq("vertices.csv"));
        v << "no,x,y\n";
        for (int64_t i = 0; i < g->NNodes(); i++)
            v << i << "," << g->NodeVec()[i].Coord(0) << "," << g->NodeVec()[i].Coord(1) << "\n";
        std::ofstream t(Arq("triangles.csv")), p(Arq("pontos.csv"));
        t << "n0,n1,n2\n";
        p << "tri,x,y,peso\n";
        int64_t linha = 0;
        ForEachPoint(fModel, [&](TPZCompEl *cel, const TPZVec<REAL> &x, REAL peso, int ip, int64_t) {
            if (ip == 0) { // primeiro ponto do elemento: uma linha nova em triangles.csv
                TPZGeoEl *gel = cel->Reference();
                t << gel->NodeIndex(0) << "," << gel->NodeIndex(1) << "," << gel->NodeIndex(2) << "\n";
                fElementos.push_back(cel);
                linha++;
            }
            p << linha - 1 << "," << x[0] << "," << x[1] << "," << peso << "\n";
        });
        for (const char *n : {"snaps_pw.csv", "snaps_pa.csv", "snaps_sw.csv"}) std::ofstream(Arq(n)); // vazios
        std::ofstream(Arq("fs.csv")) << "t_h,FS\n";
        std::ofstream(Arq("srm_nos.csv")) << "t_h,no,ux,uy,ux_geo,uy_geo\n";
        std::ofstream(Arq("srm_pontos.csv")) << "t_h,ponto,raizJ2_epsp,hardening,tipo\n";
        std::ofstream(Arq("nivel.csv")) << "t_h,nivel\n";
    }

    // um instante do transiente: p_w por no, o nivel do reservatorio e a saturacao (pela curva de retencao, p_a = 0)
    void Snapshot(REAL t, const std::vector<REAL> &pw, REAL nivel, const TPZUnsaturatedHydraulics &h)
    {
        Snapshot(t, pw, h);
        std::ofstream(Arq("nivel.csv"), std::ios::app) << t << "," << nivel << "\n";
    }
    // o mesmo sem reservatorio (chuva)
    void Snapshot(REAL t, const std::vector<REAL> &pw, const TPZUnsaturatedHydraulics &h)
    {
        std::ofstream fw(Arq("snaps_pw.csv"), std::ios::app), fa(Arq("snaps_pa.csv"), std::ios::app),
            fs(Arq("snaps_sw.csv"), std::ios::app);
        fw << t, fa << t, fs << t;
        for (REAL p : pw) fw << "," << p, fa << ",0", fs << "," << h.Evaluate(p, 0.).S;
        fw << "\n", fa << "\n", fs << "\n";
    }

    // o fator de seguranca do instante e o estado do fim do GIM (o que esta carregado na malha e na memoria)
    void FactorOfSafety(REAL t, REAL fs)
    {
        std::ofstream(Arq("fs.csv"), std::ios::app) << t << "," << fs << "\n";
        TMat *mat = fModel.Material();
        const int iu = mat->VariableIndex("Displacement");
        std::vector<REAL> ux(fModel.GeoMesh()->NNodes(), 0.), uy(ux);
        for (TPZCompEl *cel : fElementos) {
            TPZGeoEl *gel = cel->Reference();
            for (int i = 0; i < 3; i++) {
                TPZManVector<REAL, 3> qsi(2, 0.);
                gel->CenterPoint(i, qsi);
                TPZManVector<STATE, 9> u;
                cel->Solution(qsi, iu, u);
                ux[gel->NodeIndex(i)] = u[0], uy[gel->NodeIndex(i)] = u[1];
            }
        }
        std::ofstream n(Arq("srm_nos.csv"), std::ios::app);
        for (size_t i = 0; i < ux.size(); i++) n << t << "," << i << "," << ux[i] << "," << uy[i] << ",0,0\n";
        std::ofstream p(Arq("srm_pontos.csv"), std::ios::app);
        int64_t k = 0;
        ForEachPoint(fModel, [&](TPZCompEl *, const TPZVec<REAL> &, REAL, int, int64_t im) {
            PlasticState(p, t, k++, mat->MemItem(im));
        });
    }

    // O fator de seguranca do instante com o estado do fim do SRM de outro modelo, numa malha diferente desta (a do
    // fator, refinada de novo a cada instante; an: a analise dele): o deslocamento em cada vertice desta malha sai do
    // elemento do fator que contem o vertice, no fim do caminho e no estado geostatico (an.WeightState(): o
    // animacao.py desenha a diferenca, o mecanismo); o estado plastico de cada ponto de integracao desta malha, do ponto
    // de integracao do fator mais proximo, no elemento que o contem
    template <class TAnalysis>
    void FactorOfSafety(REAL t, REAL fs, UnsatModel<TMat> &fator, TAnalysis &an)
    {
        std::ofstream(Arq("fs.csv"), std::ios::app) << t << "," << fs << "\n";
        TMat *mat = fator.Material();
        const PointLocator loc(fator.GeoMesh());
        // elemento comp e pontos de integracao (x, y, indice da memoria) de cada elemento do fator, pelo indice geometrico
        std::map<int64_t, TPZCompEl *> cels;
        std::map<int64_t, std::vector<std::array<REAL, 3>>> pontos;
        ForEachPoint(fator, [&](TPZCompEl *cel, const TPZVec<REAL> &x, REAL, int, int64_t im) {
            cels[cel->Reference()->Index()] = cel;
            pontos[cel->Reference()->Index()].push_back({x[0], x[1], REAL(im)});
        });
        const int iu = mat->VariableIndex("Displacement");
        TPZGeoMesh *g = fModel.GeoMesh();
        // deslocamento nos vertices desta malha com a solucao que estiver carregada na malha do fator
        auto deslocamentos = [&]() {
            std::vector<std::array<STATE, 2>> d(g->NNodes(), {0., 0.});
            for (int64_t i = 0; i < g->NNodes(); i++) {
                TPZManVector<REAL, 3> qsi(2, 0.);
                TPZGeoEl *gel = loc.Find(g->NodeVec()[i].Coord(0), g->NodeVec()[i].Coord(1), qsi);
                TPZManVector<STATE, 9> v(2, 0.);
                if (gel && cels.count(gel->Index())) cels[gel->Index()]->Solution(qsi, iu, v);
                d[i] = {v[0], v[1]};
            }
            return d;
        };
        const auto u = deslocamentos(); // fim do caminho (carregado pelo SRMArcLength ou pelo GIMArcLength)
        std::vector<std::array<STATE, 2>> ugeo(g->NNodes(), {0., 0.}); // GIM: sem estado geostatico, u total
        if (an.WeightState().Rows() == an.PathState().Rows()) {
            an.LoadSolution(an.WeightState());
            ugeo = deslocamentos();
            an.LoadSolution(an.PathState()); // a malha do fator volta ao fim do caminho
        }
        std::ofstream n(Arq("srm_nos.csv"), std::ios::app);
        for (int64_t i = 0; i < g->NNodes(); i++)
            n << t << "," << i << "," << u[i][0] << "," << u[i][1] << "," << ugeo[i][0] << "," << ugeo[i][1] << "\n";
        std::ofstream p(Arq("srm_pontos.csv"), std::ios::app);
        int64_t k = 0;
        ForEachPoint(fModel, [&](TPZCompEl *, const TPZVec<REAL> &x, REAL, int, int64_t) {
            TPZManVector<REAL, 3> qsi(2, 0.);
            TPZGeoEl *gel = loc.Find(x[0], x[1], qsi);
            if (!gel || !pontos.count(gel->Index())) DebugStop(); // ponto desta malha fora da malha do fator
            const std::array<REAL, 3> *perto = nullptr;
            REAL d2 = 1.e300;
            for (const auto &q : pontos[gel->Index()]) {
                const REAL d = (q[0] - x[0]) * (q[0] - x[0]) + (q[1] - x[1]) * (q[1] - x[1]);
                if (d < d2) d2 = d, perto = &q;
            }
            PlasticState(p, t, k++, mat->MemItem(int64_t((*perto)[2])));
        });
    }

    void Summary(REAL tRebaixamento, REAL nivelInicial, REAL nivelFinal)
    {
        std::ofstream s(Arq("summary.txt"));
        s << "trebaix " << tRebaixamento << "\nnivel0 " << nivelInicial << "\nnivel1 " << nivelFinal << "\n";
    }
    // chuva: intensidade (mm/h) e duracao (h)
    void RainSummary(REAL intensidade, REAL tChuva)
    {
        std::ofstream(Arq("summary.txt")) << "I " << intensidade << "\ntrain " << tChuva << "\n";
    }

private:
    std::string Arq(const std::string &n) const { return fPasta + "/" + n; }

    // uma linha de srm_pontos.csv: sqrt(J2(eps_p)), hardening e tipo de resposta do ponto
    template <class TMem>
    static void PlasticState(std::ofstream &p, REAL t, int64_t k, const TMem &mem)
    {
        const auto &e = mem.m_elastoplastic_state;
        TPZTensor<REAL> eps(e.m_eps_p); // deformacoes de corte de engenharia -> componentes do tensor
        eps.XY() *= 0.5, eps.XZ() *= 0.5, eps.YZ() *= 0.5;
        p << t << "," << k << "," << std::sqrt(std::max(REAL(0), eps.J2())) << "," << e.m_hardening << ","
          << e.m_m_type << "\n";
    }

public:
    // cada ponto de integracao dos elementos do dominio do modelo, na ordem da memoria: f(cel, x, peso, ip, indice do
    // item de memoria do ponto). Publica: o montecarlo.cpp usa para descrever o mecanismo
    template <class F>
    static void ForEachPoint(UnsatModel<TMat> &model, F f)
    {
        TMat *mat = model.Material();
        for (TPZCompEl *cel : model.CompMesh()->ElementVec()) {
            if (!cel || !cel->Reference() || cel->Reference()->MaterialId() != mat->Id()) continue;
            TPZGeoEl *gel = cel->Reference();
            if (gel->NCornerNodes() != 3) DebugStop(); // o animacao.py so desenha triangulos
            auto *mf = dynamic_cast<TPZMultiphysicsElement *>(cel);
            TPZManVector<int, 3> ordens(mf->NMeshes(), 0);
            for (int i = 0; i < mf->NMeshes(); i++)
                if (auto *sp = dynamic_cast<TPZInterpolationSpace *>(mf->Element(i))) ordens[i] = sp->MaxOrder();
            const int o = mat->IntegrationRuleOrder(ordens);
            std::unique_ptr<TPZIntPoints> regra(gel->CreateSideIntegrationRule(gel->NSides() - 1, o));
            TPZManVector<int, 3> ordem(2, o);
            regra->SetOrder(ordem);
            TPZManVector<int64_t> mem;
            cel->GetMemoryIndices(mem);
            for (int ip = 0; ip < regra->NPoints(); ip++) {
                TPZManVector<REAL, 3> qsi(2), x(3), xsi(2);
                REAL w;
                regra->Point(ip, qsi, w);
                gel->X(qsi, x);
                TPZFNMatrix<9, REAL> jac(2, 2), axes(2, 3), jacinv(2, 2);
                REAL detjac;
                gel->Jacobian(qsi, jac, axes, detjac, jacinv);
                f(cel, x, w * std::abs(detjac), ip, mem[ip]);
            }
        }
    }

private:
    std::string fPasta;
    UnsatModel<TMat> &fModel;
    std::vector<TPZCompEl *> fElementos; // os elementos na ordem de triangles.csv
};
