// KLRandomField.cpp — ver KLRandomField.h

#include "KLRandomField.h"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <random>
#include <set>
#include <stdexcept>

#include "TPZBFileStream.h"
#include "TPZLapackEigenSolver.h"
#include "TPZLinearAnalysis.h"
#include "TPZMatKLKernel.h"
#include "TPZSimpleTimer.h"
#include "pzdoublestrmatriz.h"
#include "pzgeoel.h"
#include "pzinterpolationspace.h"
#include "pzsbndmat.h"
#include "pzstack.h"

namespace {
const int kMatKL = 1;  // ESoil
const int kCacheVersion = 2;
}

TPZKLRandomField::TPZKLRandomField(TPZGeoMesh *gmesh, const TOptions &opt) : fOpt(opt), fGMesh(gmesh) {
    BuildMesh();
}

TPZKLRandomField::~TPZKLRandomField() {
    delete fCMesh;
}

int TPZKLRandomField::NEquations() const {
    return (int)fCMesh->NEquations();
}

void TPZKLRandomField::BuildMesh() {
    fCMesh = new TPZCompMesh(fGMesh.get());
    fCMesh->SetDimModel(2);
    fCMesh->SetDefaultOrder(fOpt.porder);
    fCMesh->SetAllCreateFunctionsContinuous();
    auto *mat = new TPZMatKLKernel(kMatKL, 2, fOpt.Lx, fOpt.Ly);
    fCMesh->InsertMaterialObject(mat);
    std::set<int> matids = {kMatKL};
    fCMesh->AutoBuild(matids);
    fCMesh->LoadReferences();
}

REAL TPZKLRandomField::VarianceError(int M) const {
    REAL s = 0.;
    for (int k = 0; k < M && k < (int)fLambda.size(); k++) s += std::max(fLambda[k], REAL(0.));
    return 1. - s / fArea;
}

bool TPZKLRandomField::ReadCache() {
    if (fOpt.cacheFile.empty()) return false;
    std::ifstream test(fOpt.cacheFile, std::ios::binary);
    if (!test.good()) return false;
    test.close();
    TPZBFileStream in;
    in.OpenRead(fOpt.cacheFile);
    int version = 0, neq = 0, porder = 0, M = 0, nl = 0;
    REAL Lx = 0., Ly = 0., area = 0.;
    in.Read(&version, 1);
    in.Read(&neq, 1);
    in.Read(&porder, 1);
    in.Read(&Lx, 1);
    in.Read(&Ly, 1);
    in.Read(&area, 1);
    if (version != kCacheVersion || neq != NEquations() || porder != fOpt.porder ||
        std::fabs(Lx - fOpt.Lx) > 1.e-12 || std::fabs(Ly - fOpt.Ly) > 1.e-12) {
        std::cout << "[KL] cache " << fOpt.cacheFile << " incompatível; recalculando\n";
        return false;
    }
    in.Read(&nl, 1);
    fLambda.resize(nl);
    if (nl) in.Read(fLambda.data(), nl);
    in.Read(&M, 1);
    fPhi.Read(in, nullptr);
    fArea = area;
    fM = M;
    if (fPhi.Rows() != neq || fPhi.Cols() != M) return false;
    std::cout << "[KL] lido de " << fOpt.cacheFile << ": " << M << " modos\n";
    return true;
}

void TPZKLRandomField::WriteCache() const {
    if (fOpt.cacheFile.empty()) return;
    TPZBFileStream out;
    out.OpenWrite(fOpt.cacheFile);
    int version = kCacheVersion, neq = NEquations(), porder = fOpt.porder, nl = (int)fLambda.size(), M = fM;
    REAL Lx = fOpt.Lx, Ly = fOpt.Ly, area = fArea;
    out.Write(&version, 1);
    out.Write(&neq, 1);
    out.Write(&porder, 1);
    out.Write(&Lx, 1);
    out.Write(&Ly, 1);
    out.Write(&area, 1);
    out.Write(&nl, 1);
    if (nl) out.Write(fLambda.data(), nl);
    out.Write(&M, 1);
    fPhi.Write(out, 0);
}

void TPZKLRandomField::Compute() {
    const int64_t n = fCMesh->NEquations();
    if (ReadCache()) return;

    TPZSimpleTimer timer("KL");
    auto *mat = dynamic_cast<TPZMatKLKernel *>(fCMesh->FindMaterial(kMatKL));
    pzdoublestrmatriz<STATE> sm(fCMesh);
    sm.SetCAssembly(pzdoublestrmatriz<STATE>::ECAssembly::Galerkin);
    TPZFMatrix<STATE> C, B, rhs;
    TPZAutoPointer<TPZGuiInterface> gui;
    mat->SetMatrixA();
    sm.Assemble(C, rhs, gui);
    mat->SetMatrixB();
    sm.Assemble(B, rhs, gui);
    mat->SetMatrixA();
    // |Ω| = Σ_e Σ_ip w detJ (as funções H1 hierárquicas do NeoPZ não formam partição da unidade: Σ B_ij ≠ |Ω|)
    fArea = 0.;
    for (TPZCompEl *cel : fCMesh->ElementVec()) {
        auto *intel = dynamic_cast<TPZInterpolationSpace *>(cel);
        if (!intel || intel->Dimension() != 2) continue;
        const TPZIntPoints &rule = intel->GetIntegrationRule();
        TPZMaterialDataT<STATE> data;
        intel->InitMaterialData(data);
        TPZManVector<REAL, 3> qsi(2, 0.);
        REAL w;
        for (int ip = 0; ip < rule.NPoints(); ip++) {
            rule.Point(ip, qsi, w);
            intel->ComputeRequiredData(data, qsi);
            fArea += w * data.detjac;
        }
    }
    std::cout << "[KL] " << n << " equações, |Ω| = " << fArea << ", C e B montadas em " << timer.ReturnTimeDouble() / 1000.
              << " s\n";

    // problema simétrico-definido em armazenamento de banda cheia (dsbgv)
    TPZSBMatrix<STATE> Csb(n, n - 1), Bsb(n, n - 1);
    for (int64_t i = 0; i < n; i++)
        for (int64_t j = i; j < n; j++) {
            Csb.PutVal(i, j, 0.5 * (C(i, j) + C(j, i)));
            Bsb.PutVal(i, j, 0.5 * (B(i, j) + B(j, i)));
        }
    C.Resize(0, 0);
    TPZLapackEigenSolver<STATE> solver;
    TPZVec<CSTATE> w;
    TPZFMatrix<CSTATE> V;
    const int info = solver.SolveGeneralisedEigenProblem(Csb, Bsb, w, V);
    if (info != 0) throw std::runtime_error("KL: dsbgv falhou");
    // dsbgv: autovalores crescentes -> ordem decrescente
    fLambda.resize(n);
    for (int64_t k = 0; k < n; k++) fLambda[k] = std::real(w[n - 1 - k]);

    if (fOpt.targetVarianceError > 0.) {
        fM = (int)n;
        for (int M = 1; M <= n; M++)
            if (VarianceError(M) <= fOpt.targetVarianceError) {
                fM = M;
                break;
            }
    } else {
        fM = (fOpt.nModes > 0) ? (int)std::min<int64_t>(fOpt.nModes, n) : (int)n;
    }
    // Φ_k = √λ_k φ_k, com φ_k normalizado em L2 (φᵀ B φ = 1)
    fPhi.Redim(n, fM);
    TPZFMatrix<STATE> v(n, 1), Bv;
    for (int k = 0; k < fM; k++) {
        const int64_t col = n - 1 - k;
        for (int64_t i = 0; i < n; i++) v(i, 0) = std::real(V(i, col));
        B.Multiply(v, Bv);
        REAL nrm2 = 0.;
        for (int64_t i = 0; i < n; i++) nrm2 += v(i, 0) * Bv(i, 0);
        const REAL scale = std::sqrt(std::max(fLambda[k], REAL(0.)) / nrm2);
        for (int64_t i = 0; i < n; i++) fPhi(i, k) = v(i, 0) * scale;
    }
    std::cout << "[KL] autoproblema resolvido em " << timer.ReturnTimeDouble() / 1000. << " s; λ1 = " << fLambda[0]
              << ", M = " << fM << ", ε_M = " << VarianceError(fM) << " (todos os modos: " << VarianceError((int)n)
              << ")\n";
    WriteCache();
}

void TPZKLRandomField::LocatePoint(const TPZManVector<REAL, 3> &x, TPZInterpolationSpace *&cel,
                                   TPZManVector<REAL, 3> &qsi, int64_t &start) const {
    TPZManVector<REAL, 3> xx(x), q(2, 0.);
    if (start < 0 || start >= fGMesh->NElements() || !fGMesh->Element(start)) start = 0;
    while (start < fGMesh->NElements() &&
           (!fGMesh->Element(start) || fGMesh->Element(start)->Dimension() != 2))
        start++;
    TPZGeoEl *gel = fGMesh->FindElement(xx, q, start, 2);
    auto check = [&](TPZGeoEl *g, TPZManVector<REAL, 3> &qq) {
        if (!g || g->Dimension() != 2 || !g->Reference()) return false;
        TPZManVector<REAL, 3> y(3, 0.);
        g->X(qq, y);
        const REAL d = std::hypot(y[0] - x[0], y[1] - x[1]);
        return d < 1.e-8 && g->IsInParametricDomain(qq, 1.e-6);
    };
    if (!check(gel, q)) {
        gel = nullptr;
        for (int64_t i = 0; i < fGMesh->NElements(); i++) {
            TPZGeoEl *g = fGMesh->Element(i);
            if (!g || g->Dimension() != 2 || !g->Reference() || g->HasSubElement()) continue;
            TPZManVector<REAL, 3> qq(2, 0.);
            g->ComputeXInverse(xx, qq, 1.e-10);
            if (check(g, qq)) {
                gel = g;
                q = qq;
                start = i;
                break;
            }
        }
    }
    if (!gel) throw std::runtime_error("KL: ponto fora da malha do campo aleatório");
    cel = dynamic_cast<TPZInterpolationSpace *>(gel->Reference());
    qsi = q;
}

int TPZKLRandomField::AddTargetSet(const std::vector<TPZManVector<REAL, 3>> &points) {
    TTargetSet set;
    set.pts = points;
    const size_t np = points.size();
    set.cel.resize(np);
    set.qsi.resize(np);
    set.var.assign(np, 1.);
    int64_t start = 0;
    for (size_t i = 0; i < np; i++) LocatePoint(points[i], set.cel[i], set.qsi[i], start);
    // v(x) = Σ_k Φ_k(x)²: malha carregada com as M colunas de Φ
    fCMesh->LoadSolution(fPhi);
    TPZMaterialDataT<STATE> data;
    for (size_t i = 0; i < np; i++) {
        set.cel[i]->ComputeSolution(set.qsi[i], data, false);
        REAL v = 0.;
        for (int k = 0; k < fM; k++) v += data.sol[k][0] * data.sol[k][0];
        set.var[i] = v;
    }
    fSets.push_back(std::move(set));
    return (int)fSets.size() - 1;
}

void TPZKLRandomField::Xi(uint64_t seed, int64_t sample, int field, std::vector<REAL> &xi) const {
    std::seed_seq seq{(uint32_t)(seed & 0xffffffffu), (uint32_t)(seed >> 32), (uint32_t)(sample & 0xffffffffu),
                      (uint32_t)((uint64_t)sample >> 32), (uint32_t)field};
    std::mt19937_64 gen(seq);
    std::normal_distribution<double> N01(0., 1.);
    xi.resize(fM);
    for (int k = 0; k < fM; k++) xi[k] = (REAL)N01(gen);
}

void TPZKLRandomField::Evaluate(const std::vector<std::vector<REAL>> &xi, const std::vector<int> &sets,
                                std::vector<std::vector<std::vector<REAL>>> &values) {
    const int nf = (int)xi.size();
    TPZFMatrix<STATE> X(fM, nf), Z;
    for (int f = 0; f < nf; f++)
        for (int k = 0; k < fM; k++) X(k, f) = xi[f][k];
    fPhi.Multiply(X, Z);  // n x nf: coeficientes das realizações H_f
    fCMesh->LoadSolution(Z);
    values.assign(nf, std::vector<std::vector<REAL>>(fSets.size()));
    TPZMaterialDataT<STATE> data;
    for (int s : sets) {
        const TTargetSet &set = fSets[s];
        for (int f = 0; f < nf; f++) values[f][s].resize(set.pts.size());
        for (size_t i = 0; i < set.pts.size(); i++) {
            TPZManVector<REAL, 3> q(set.qsi[i]);
            set.cel[i]->ComputeSolution(q, data, false);
            const REAL scale = fOpt.normalizeVariance ? 1. / std::sqrt(std::max(set.var[i], REAL(1.e-12))) : 1.;
            for (int f = 0; f < nf; f++) values[f][s][i] = data.sol[f][0] * scale;
        }
    }
}

REAL TPZKLRandomField::Lognormal(REAL h, REAL mean, REAL cov) {
    if (cov <= 0.) return mean;
    const REAL s2 = std::log(1. + cov * cov);
    const REAL mu = std::log(mean) - 0.5 * s2;
    return std::exp(mu + std::sqrt(s2) * h);
}

void TPZKLRandomField::WriteVTK(const std::vector<REAL> &xi, const std::string &file, int step) {
    TPZFMatrix<STATE> X(fM, 1), Z;
    for (int k = 0; k < fM; k++) X(k, 0) = xi[k];
    fPhi.Multiply(X, Z);
    TPZLinearAnalysis an(fCMesh, false);
    an.LoadSolution(Z);
    TPZStack<std::string> scal, vec;
    scal.Push("KLField");
    an.DefineGraphMesh(2, scal, vec, file);
    an.SetStep(step);
    an.PostProcess(1);
}
