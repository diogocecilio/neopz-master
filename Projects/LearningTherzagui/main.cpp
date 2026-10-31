// LearningTherzagui: coluna de Terzaghi resolvida passo a passo (projeto de estudo).
// Por enquanto so os includes; o programa sera escrito aqui.

// geometria
#include "TPZGenGrid2D.h"
#include "TPZGeoLinear.h"
#include "pzgeoel.h"
#include "pzgeoquad.h"
#include "pzgeotriangle.h"
#include "pzgmesh.h"
#include "tpzgeoelrefpattern.h"

// malhas computacionais
#include "TPZBndCondT.h"
#include "TPZMultiphysicsCompMesh.h"
#include "TPZNullMaterial.h"
#include "pzcmesh.h"

// material u - p_w - p_a e Mohr-Coulomb PV2 do artigo
#include "Plasticity/TPZElasticResponse.h"
#include "Plasticity/TPZMatUnsaturatedUPwPa.h"
#include "Plasticity/TPZPlasticStepVoigt.h"
#include "Plasticity/TPZUnsaturatedHydraulics.h"
#include "Plasticity/TPZYCMohrCoulombPV2.h"

// analise e solver
#include "TPZLinearAnalysis.h"
#include "TPZSkylineNSymStructMatrix.h"
#include "TPZSpStructMatrix.h"
#include "pzstepsolver.h"

#include <cmath>
#include <fstream>
#include <iostream>

#include "tpzgeoblend.h"
#include "tpzarc3d.h"
#include "TPZGeoCube.h"
#include <array>
#include <vector>
#include "TPZVTKGeoMesh.h"
#include "Plasticity/TPZMatPoroElastoPlasticUP.h"
#include "pzgeoelbc.h"
#include <set>
#include "Plasticity/TPZPoroElastoPlasticUPAnalysis.h"
#include "TPZSkylineNSymStructMatrix.h"
#include "pzstepsolver.h"
using TMC = TPZPlasticStepVoigt<TPZYCMohrCoulombPV2, TPZElasticResponse>;
using TMatUP = TPZMatPoroElastoPlasticUP<TMC, TPZElastoPlasticMem>;



// ids de material da coluna de Terzaghi (malha geometrica e computacional)
enum EIds {
    matId = 1,      // volume
    bcTop = -1,     // topo (carga)
    bcBottom = -2,  // base
    bcXZ = -3,      // plano xz (y = 0)
    bcYZ = -4,      // plano yz (x = 0)
    bcOuter = -5,   // contorno circular
    matArc = -6,    // arcos TPZArc3D (auxiliar, nao registrar como material)
    bcTopP = -11    // gemeo do topo (drenagem p = 0)
};

// Parametros do problema de Terzaghi (Wang, Theory of Linear Poroelasticity, sec. 6.5)
// Unidades kN, m, s, kPa. z para baixo, topo z = 0 carregado e drenado, base z = L fixa e impermeavel.
struct PoroElasticParams {
    REAL E = 1.e4;          // modulo de Young drenado (kPa)
    REAL nu = 0.3;          // Poisson drenado
    REAL kh = 1.e-7;        // condutividade hidraulica (m/s)
    REAL gammaw = 9.81;     // peso especifico da agua (kN/m3)
    REAL alpha = 1.;        // coeficiente de Biot
    REAL invM = 0.;         // 1/M de Biot (kPa^-1); 0: constituintes incompressiveis
    REAL q = 10.;           // sigma_o: carga no topo (kPa, compressao)
    REAL L = 1.;            // comprimento da coluna (m)
    REAL phi = 30. * M_PI / 180.; // so para o Mohr-Coulomb
    REAL c = 1.e10;         // coesao enorme: solo elastico (nao confundir com o c de Wang = Cv())

    // mobilidade k = kh / gammaw (m^2/(kPa s))
    REAL Mobility() const { return kh / gammaw; }
    // K_v: modulo edometrico drenado
    REAL Kv() const { return E * (1. - nu) / ((1. + nu) * (1. - 2. * nu)); }
    // K_v^(u) = K_v + alpha^2 M (infinito se invM = 0)
    REAL KvU() const { return invM > 0. ? Kv() + alpha * alpha / invM : INFINITY; }
    // gamma: eficiencia de carregamento uniaxial, p0 = gamma sigma_o
    REAL Gamma() const { return alpha / (Kv() * invM + alpha * alpha); }
    // S: armazenamento uniaxial
    REAL Storage() const { return invM + alpha * alpha / Kv(); }
    // c de Wang: coeficiente de adensamento c = k / S (m^2/s)
    REAL Cv() const { return Mobility() / Storage(); }
    // tau = c t / L^2 e o tempo correspondente
    REAL Tau(REAL t) const { return Cv() * t / (L * L); }
    REAL Time(REAL tau) const { return tau * L * L / Cv(); }

    // p(z,t), eq. (6.34)
    REAL Pressure(REAL z, REAL t) const {
        const REAL p0 = Gamma() * q;
        if (t <= 0.) return p0;
        REAL s = 0.;
        for (int m = 0; m < 2000; m++) {
            const REAL a = (2 * m + 1) * M_PI / (2. * L);
            const REAL e = std::exp(-a * a * Cv() * t);
            if (e < 1.e-16) break;
            s += e * std::sin(a * z) / (2 * m + 1);
        }
        return 4. * p0 / M_PI * s;
    }

    // w(z,t) para baixo: dw/dz = -(sigma_o - alpha p)/K_v com w(L) = 0.
    // t = 0+ da (6.37) sigma_o (L - z)/K_v^(u); t -> infinito da sigma_o (L - z)/K_v
    REAL Settlement(REAL z, REAL t) const {
        const REAL p0 = Gamma() * q;
        REAL s = 0.;
        for (int m = 0; m < 2000; m++) {
            const REAL a = (2 * m + 1) * M_PI / (2. * L);
            const REAL e = (t <= 0.) ? 1. : std::exp(-a * a * Cv() * t);
            if (e < 1.e-16) break;
            s += e * std::cos(a * z) / ((2 * m + 1) * a);
        }
        const REAL intP = 4. * p0 / M_PI * s; // integral de z a L de p
        return (q * (L - z) - alpha * intP) / Kv();
    }
};


// Material u-p com solucao exata para pos-processamento (como o caso ExactPressure do TPZDarcyFlow).
// A solucao exata depende do tempo: fExact(x, t, u, p), com t = fTime atualizado a cada passo.
class TMatUPExact : public TMatUP {
public:
    using ExactUP = std::function<void(const TPZVec<REAL> &x, REAL t, TPZVec<STATE> &u, STATE &p)>;
    // indices >= 100 sao reservados pelo NeoPZ (solucao por elemento)
    enum EExactVar { EExactPorePressure = 50, EExactDisplacement, EErrorPorePressure, EErrorDisplacement };

    TMatUPExact(int id, EKinematics kin) : TMatUP(id, kin) {}
    TMatUPExact(const TMatUPExact &cp) = default;

    void SetExactSolution(ExactUP f) { fExact = f; }
    void SetTime(REAL t) { fTime = t; }
    REAL Time() const { return fTime; }

    TPZMaterial *NewMaterial() const override { return new TMatUPExact(*this); }

    int VariableIndex(const std::string &name) const override {
        if (name == "ExactPorePressure") return EExactPorePressure;
        if (name == "ExactDisplacement") return EExactDisplacement;
        if (name == "ErrorPorePressure") return EErrorPorePressure;
        if (name == "ErrorDisplacement") return EErrorDisplacement;
        return TMatUP::VariableIndex(name);
    }

    int NSolutionVariables(int var) const override {
        if (var == EExactPorePressure || var == EErrorPorePressure) return 1;
        if (var == EExactDisplacement || var == EErrorDisplacement) return 3;
        return TMatUP::NSolutionVariables(var);
    }

    void Solution(const TPZVec<TPZMaterialDataT<STATE>> &datavec, int var, TPZVec<STATE> &sol) override {
        if (var < EExactPorePressure || var > EErrorDisplacement) {
            TMatUP::Solution(datavec, var, sol);
            return;
        }
        if (!fExact) DebugStop(); // solucao exata nao informada
        TPZManVector<STATE, 3> ue(3, 0.);
        STATE pe = 0.;
        fExact(datavec[0].x, fTime, ue, pe);
        sol.Resize(NSolutionVariables(var));
        switch (var) {
            case EExactPorePressure: sol[0] = pe; break;
            case EExactDisplacement: for (int i = 0; i < 3; i++) sol[i] = ue[i]; break;
            case EErrorPorePressure: sol[0] = datavec[1].sol[0][0] - pe; break;
            case EErrorDisplacement: for (int i = 0; i < 3; i++) sol[i] = datavec[0].sol[0][i] - ue[i]; break;
        }
    }

private:
    ExactUP fExact;
    REAL fTime = 0.;
};



TPZGeoMesh *CreateQuarterCylinderMesh(REAL R, REAL L, int nCirc, int nRad,
                                      int nz, int matId, int matArc);

TPZCompMesh *CreateAtomicMesh(TPZGeoMesh *gmesh, int nstate, int order, const std::set<int> &bcids);


TPZMultiphysicsCompMesh *CreateTerzaghiCompMesh(TPZGeoMesh *gmesh, int porder, const PoroElasticParams &prm);

TPZPoroElastoPlasticUPAnalysis *CreateAnalysis(TPZMultiphysicsCompMesh *cmesh, int nthreads);

void PostProcessVTK(TPZPoroElastoPlasticUPAnalysis *an, TPZMultiphysicsCompMesh *cmesh, int step, REAL time,
                    const std::string &file, int resolution);
void WriteAxisCSV(TPZMultiphysicsCompMesh *cmesh, const PoroElasticParams &prm, int step, REAL time, std::ofstream &csv);
int main() {


    PoroElasticParams prm;
    TPZGeoMesh *gmesh = CreateQuarterCylinderMesh(0.2, prm.L, 2, 2, 20, 1, -6);
    TPZMultiphysicsCompMesh *cmesh = CreateTerzaghiCompMesh(gmesh, 2, prm);
    std::ofstream out("QuarterCylinder.vtk");
    TPZVTKGeoMesh::PrintGMeshVTK(gmesh, out, true);
    TPZPoroElastoPlasticUPAnalysis *an = CreateAnalysis(cmesh, 4);


    std::ofstream csv("Terzaghi_eixo.csv");
    csv << "step,t,tau,z,p,p_exact,w,w_exact\n";

    auto monitor = [&](int k, const TPZPoroElastoPlasticUPAnalysis::TLoadState &s) {
        std::cout << "passo " << k << "  t = " << s.fTime << " s  tau = " << prm.Cv() * s.fTime / (prm.L * prm.L)
        << std::endl;
        PostProcessVTK(an, cmesh, k, s.fTime, "Terzaghi.vtk", 1);
        WriteAxisCSV(cmesh, prm, k, s.fTime, csv);
    };
    std::vector<TPZPoroElastoPlasticUPAnalysis::TLoadState> steps;
    steps.push_back({0., 1., 0.}); // t = 0: carga aplicada, passo nao drenado
    for (REAL tau : {0.01, 0.04, 0.1, 0.4, 1.0}) {
        const REAL t0 = steps.back().fTime, t1 = prm.Time(tau);
        const int n = 10; // passos de tempo em cada trecho
        for (int i = 1; i <= n; i++) steps.push_back({t0 + (t1 - t0) * i / n, 1., 0.});
    }
    const bool ok = an->Run(steps, monitor);
    if (!ok) std::cout << "Run falhou: um passo nao convergiu apos as bissecoes" << std::endl;
    return 0;
}

// Pos-processamento VTK de um passo: atualiza o tempo da solucao exata e grava
// <file sem .vtk>.scal_vec.<step>.vtk
//   resolution: 0 = so os vertices; 1 = cada elemento subdividido (melhor para u quadratico)
void PostProcessVTK(TPZPoroElastoPlasticUPAnalysis *an, TPZMultiphysicsCompMesh *cmesh, int step, REAL time,
                    const std::string &file, int resolution)
{
    auto *mat = dynamic_cast<TMatUPExact *>(cmesh->FindMaterial(matId));
    if (mat) mat->SetTime(time);
    TPZManVector<std::string, 4> scal = {"PorePressure", "ExactPorePressure", "ErrorPorePressure", "DisplacementZ"};
    TPZManVector<std::string, 3> vec = {"Displacement", "ExactDisplacement", "ErrorDisplacement"};
    an->DefineGraphMesh(3, scal, vec, file);
    an->SetStep(step);
    an->PostProcess(resolution, 3);
}

// Analise incremental u-p: matriz skyline nao simetrica + LU (sem pivotamento).
// O driver numera as equacoes de pressao depois das de deslocamento (sem otimizar banda),
// o que mantem a LU estavel no passo nao drenado com constituintes incompressiveis.
TPZPoroElastoPlasticUPAnalysis *CreateAnalysis(TPZMultiphysicsCompMesh *cmesh, int nthreads)
{
    auto *mat = dynamic_cast<TPZMatPoroElastoPlasticUPBase *>(cmesh->FindMaterial(matId));
    if (!mat) DebugStop();
    auto *an = new TPZPoroElastoPlasticUPAnalysis(cmesh, mat);

    // matriz estrutural
    TPZSkylineNSymStructMatrix<STATE> skyl(cmesh);
    skyl.SetNumThreads(nthreads);
    an->SetStructuralMatrix(skyl);

    // solver direto LU
    TPZStepSolver<STATE> step;
    step.SetDirect(ELU);
    an->SetSolver(step);

    // Newton: tolerancia do residuo normalizado, max. iteracoes, max. bissecoes
    an->SetNewtonParameters(1.e-8, 10, 4);
    // norma do residuo na base hierarquica do NeoPZ (a nodal supoe Hex20 serendipity)
    an->SetNodalResidualNorm(false);
    an->SetVerbose(1); // uma linha por passo
    return an;
}

// Malha atomica H1 com material nulo: nstate = 3 (deslocamento) ou 1 (poropressao)
TPZCompMesh *CreateAtomicMesh(TPZGeoMesh *gmesh, int nstate, int order, const std::set<int> &bcids)
{
    auto *cmesh = new TPZCompMesh(gmesh);
    cmesh->SetDimModel(3);
    cmesh->SetDefaultOrder(order);
    cmesh->SetAllCreateFunctionsContinuous();
    auto *mat = new TPZNullMaterial<STATE>(1, 3, nstate);
    cmesh->InsertMaterialObject(mat);
    TPZFMatrix<STATE> val1(nstate, nstate, 0.);
    TPZManVector<STATE, 3> val2(nstate, 0.);
    for (int id : bcids) cmesh->InsertMaterialObject(mat->CreateBC(mat, id, 0, val1, val2));
    cmesh->AutoBuild();
    return cmesh;
}


TPZMultiphysicsCompMesh *CreateTerzaghiCompMesh(TPZGeoMesh *gmesh, int porder, const PoroElasticParams &prm)
{
    using B = TPZMatPoroElastoPlasticUPBase;

    // gemeo do topo: o topo precisa de uma condicao em u (carga) e outra em p (drenado)
    const int64_t nel = gmesh->NElements();
    for (int64_t i = 0; i < nel; i++) {
        TPZGeoEl *gel = gmesh->Element(i);
        if (gel && !gel->HasSubElement() && gel->MaterialId() == bcTop)
            TPZGeoElBC(gel, gel->NSides() - 1, bcTopP);
    }

    const std::set<int> bcids = {bcTop, bcBottom, bcXZ, bcYZ, bcOuter, bcTopP};
    TPZManVector<TPZCompMesh *, 2> meshvec = {CreateAtomicMesh(gmesh, 3, porder, bcids),
                                              CreateAtomicMesh(gmesh, 1, porder - 1, bcids)};

    // material de Biot, sem peso proprio; coesao enorme: o solo fica elastico
    TPZElasticResponse ER;
    ER.SetEngineeringData(prm.E, prm.nu);
    //auto *mat = new TMatUP(matId, B::EThreeDimensional);
    auto *mat = new TMatUPExact(matId, B::EThreeDimensional);
    // solucao de Wang (6.34) e recalque; prm copiado para dentro da lambda
    mat->SetExactSolution([prm](const TPZVec<REAL> &x, REAL t, TPZVec<STATE> &u, STATE &p) {
        u[0] = 0.;
        u[1] = 0.;
        u[2] = prm.Settlement(x[2], t);
        p = prm.Pressure(x[2], t);
    });
    mat->SetPlasticModel(TMC(TPZYCMohrCoulombPV2(prm.phi, prm.phi, prm.c, ER), ER));
    mat->SetBiot(prm.alpha, prm.invM);
    mat->SetPermeability(prm.Mobility());

    auto *mphys = new TPZMultiphysicsCompMesh(gmesh);
    mphys->SetDimModel(3);
    mphys->InsertMaterialObject(mat);

    TPZFMatrix<STATE> val1(3, 3, 0.);
    TPZManVector<STATE, 3> val2(3, 0.);

    // topo: tracao t = sigma n, n = (0,0,-1) e sigma_zz = -q  ->  t = (0,0,q) (empurra para baixo)
    val2[2] = prm.q;
    mphys->InsertMaterialObject(mat->CreateBC(mat, bcTop, B::ENeumannU, val1, val2));
    val2.Fill(0.);
    // topo drenado
    mphys->InsertMaterialObject(mat->CreateBC(mat, bcTopP, B::EDirichletP, val1, val2));
    // base engastada
    mphys->InsertMaterialObject(mat->CreateBC(mat, bcBottom, B::EDirichletU, val1, val2));
    // plano xz (y = 0): u_y = 0
    val1.Zero();
    val1(1, 1) = 1.;
    mphys->InsertMaterialObject(mat->CreateBC(mat, bcXZ, B::EDirichletUDirectional, val1, val2));
    // plano yz (x = 0): u_x = 0
    val1.Zero();
    val1(0, 0) = 1.;
    mphys->InsertMaterialObject(mat->CreateBC(mat, bcYZ, B::EDirichletUDirectional, val1, val2));
    // lateral: u_x = u_y = 0
    val1.Zero();
    val1(0, 0) = 1.;
    val1(1, 1) = 1.;
    mphys->InsertMaterialObject(mat->CreateBC(mat, bcOuter, B::EDirichletUDirectional, val1, val2));

    TPZManVector<int, 2> active(2, 1);
    mphys->BuildMultiphysicsSpaceWithMemory(active, meshvec, {matId}, bcids);
    return mphys;
}


TPZGeoMesh *CreateQuarterCylinderMesh(REAL R, REAL L, int nCirc, int nRad,
                                      int nz, int matId, int matArc)
{

    const REAL a = 0.5 * R;
    const REAL tol = 1.e-10 * R;

    // ---- 1. malha 2D da planta (nos unicos + quadrilateros) ----
    std::vector<std::array<REAL, 2>> xy;
    auto addNode2D = [&](REAL x, REAL y) -> int {
        for (int i = 0; i < (int)xy.size(); i++) {
            if (std::fabs(xy[i][0] - x) < tol && std::fabs(xy[i][1] - y) < tol) return i;
        }
        xy.push_back({x, y});
        return (int)xy.size() - 1;
    };
    std::vector<std::array<int, 4>> quads;
    // grade (ni+1)x(nj+1) de nos de um bloco -> quadrilateros
    auto addBlock = [&](std::vector<std::vector<int>> &g, int ni, int nj) {
        for (int i = 0; i < ni; i++) {
            for (int j = 0; j < nj; j++) {
                std::array<int, 4> q = {g[i][j], g[i + 1][j], g[i + 1][j + 1], g[i][j + 1]};
                // garante sentido anti-horario (jacobiano positivo no hexaedro)
                REAL area = 0.;
                for (int k = 0; k < 4; k++) {
                    const auto &p = xy[q[k]], &r = xy[q[(k + 1) % 4]];
                    area += p[0] * r[1] - r[0] * p[1];
                }
                if (area < 0.) std::swap(q[1], q[3]);
                quads.push_back(q);
            }
        }
    };

    // bloco central: quadrado [0,a]x[0,a]
    {
        std::vector<std::vector<int>> g(nCirc + 1, std::vector<int>(nCirc + 1));
        for (int i = 0; i <= nCirc; i++)
            for (int j = 0; j <= nCirc; j++) g[i][j] = addNode2D(a * i / nCirc, a * j / nCirc);
            addBlock(g, nCirc, nCirc);
    }
    // blocos curvos: i radial (quadrado -> arco), j angular
    for (int b = 0; b < 2; b++) {
        std::vector<std::vector<int>> g(nRad + 1, std::vector<int>(nCirc + 1));
        for (int j = 0; j <= nCirc; j++) {
            const REAL t = (REAL)j / nCirc;
            const REAL theta = M_PI / 4. * (b + t);
            // ponto no lado do quadrado: x = a (bloco 0) ou y = a (bloco 1)
            const REAL xi = (b == 0) ? a : a * (1. - t);
            const REAL yi = (b == 0) ? a * t : a;
            const REAL xo = R * std::cos(theta), yo = R * std::sin(theta);
            for (int i = 0; i <= nRad; i++) {
                const REAL s = (REAL)i / nRad;
                g[i][j] = addNode2D((1. - s) * xi + s * xo, (1. - s) * yi + s * yo);
            }
        }
        addBlock(g, nRad, nCirc);
    }

    // ---- 2. nos 3D: camadas k = 0..nz, z = L k / nz ----
    const int n2D = (int)xy.size();
    TPZGeoMesh *gmesh = new TPZGeoMesh;
    gmesh->SetDimension(3);
    auto addNode3D = [&](REAL x, REAL y, REAL z) -> int64_t {
        const int64_t id = gmesh->NodeVec().AllocateNewElement();
        TPZManVector<REAL, 3> co = {x, y, z};
        gmesh->NodeVec()[id].Initialize(co, *gmesh);
        return id;
    };
    for (int k = 0; k <= nz; k++) {
        const REAL z = L * k / nz;
        for (int i = 0; i < n2D; i++) addNode3D(xy[i][0], xy[i][1], z);
    }
    auto node3D = [&](int i2D, int k) -> int64_t { return (int64_t)k * n2D + i2D; };

    auto onY0 = [&](int i) { return std::fabs(xy[i][1]) < tol; };
    auto onX0 = [&](int i) { return std::fabs(xy[i][0]) < tol; };
    auto onArc = [&](int i) { return std::fabs(std::hypot(xy[i][0], xy[i][1]) - R) < tol; };

    using Hexa = TPZGeoElRefPattern<pzgeom::TPZGeoBlend<pzgeom::TPZGeoCube>>;
    using Quad = TPZGeoElRefPattern<pzgeom::TPZGeoBlend<pzgeom::TPZGeoQuad>>;
    using Arc = TPZGeoElRefPattern<pzgeom::TPZArc3D>;

    // ---- 3. hexaedros: nos 0-3 na camada k, nos 4-7 na camada k+1 ----
    for (int k = 0; k < nz; k++) {
        for (const auto &q : quads) {
            TPZManVector<int64_t, 8> nodes(8);
            for (int c = 0; c < 4; c++) {
                nodes[c] = node3D(q[c], k);
                nodes[c + 4] = node3D(q[c], k + 1);
            }
            new Hexa(nodes, matId, *gmesh);
        }
    }

    // ---- 4. topo (-1) e base (-2) ----
    for (const auto &q : quads) {
        TPZManVector<int64_t, 4> top(4), bottom(4);
        for (int c = 0; c < 4; c++) {
            top[c] = node3D(q[c], 0);
            bottom[c] = node3D(q[c], nz);
        }
        new Quad(top, bcTop, *gmesh);
        new Quad(bottom, bcBottom, *gmesh);
    }

    // ---- 5. faces laterais: arestas de contorno da planta x camadas ----
    for (const auto &q : quads) {
        for (int e = 0; e < 4; e++) {
            const int n0 = q[e], n1 = q[(e + 1) % 4];
            int bc = 0;
            if (onY0(n0) && onY0(n1)) bc = bcXZ;
            else if (onX0(n0) && onX0(n1)) bc = bcYZ;
            else if (onArc(n0) && onArc(n1)) bc = bcOuter;
            else continue;

            for (int k = 0; k < nz; k++) {
                TPZManVector<int64_t, 4> face = {node3D(n0, k), node3D(n1, k), node3D(n1, k + 1),
                    node3D(n0, k + 1)};
                    new Quad(face, bc, *gmesh);
            }
            if (bc != bcOuter) continue;

            // arcos TPZArc3D nas arestas horizontais do contorno circular:
            // nos (inicio, fim, ponto intermediario sobre o arco)
            const REAL t0 = std::atan2(xy[n0][1], xy[n0][0]);
            const REAL t1 = std::atan2(xy[n1][1], xy[n1][0]);
            const REAL tm = 0.5 * (t0 + t1);
            for (int k = 0; k <= nz; k++) {
                const int64_t mid = addNode3D(R * std::cos(tm), R * std::sin(tm), L * k / nz);
                TPZManVector<int64_t, 3> arc = {node3D(n0, k), node3D(n1, k), mid};
                new Arc(arc, matArc, *gmesh);
            }
        }
    }

    gmesh->BuildConnectivity();
    return gmesh;
}

// Grava no CSV os valores nodais no eixo x = y = 0 (na base hierarquica do NeoPZ o conector
// de vertice guarda o valor nodal): step, t, tau, z, p, p_exata, w, w_exata
void WriteAxisCSV(TPZMultiphysicsCompMesh *cmesh, const PoroElasticParams &prm, int step, REAL time, std::ofstream &csv)
{
    std::map<REAL, std::array<REAL, 2>> vals; // z -> (p, w)
    TPZVec<TPZCompMesh *> &mv = cmesh->MeshVector();
    for (int sp = 0; sp < 2; sp++) { // 0: deslocamento, 1: pressao
        TPZFMatrix<STATE> &sol = mv[sp]->Solution();
        for (TPZCompEl *cel : mv[sp]->ElementVec()) {
            if (!cel || !cel->Reference() || cel->Reference()->Dimension() != 3) continue;
            TPZGeoEl *gel = cel->Reference();
            for (int i = 0; i < gel->NCornerNodes(); i++) {
                TPZManVector<REAL, 3> x(3);
                gel->Node(i).GetCoordinates(x);
                if (std::hypot(x[0], x[1]) > 1.e-10) continue;
                const int64_t pos = mv[sp]->Block().Position(cel->Connect(i).SequenceNumber());
                vals[x[2]][sp == 0 ? 1 : 0] = sol(sp == 0 ? pos + 2 : pos, 0);
            }
        }
    }
    const REAL tau = prm.Cv() * time / (prm.L * prm.L);
    for (auto &[z, v] : vals)
        csv << step << ',' << time << ',' << tau << ',' << z << ',' << v[0] << ',' << prm.Pressure(z, time) << ','
        << v[1] << ',' << prm.Settlement(z, time) << '\n';
    csv.flush();
}
