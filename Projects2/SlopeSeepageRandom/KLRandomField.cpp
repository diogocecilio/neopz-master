// KLRandomField.cpp — ver KLRandomField.h

#include "KLRandomField.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <iostream>
#include <random>
#include <set>
#include <stdexcept>

#include "FileIO.h"
#include "SlopeGeometry.h"
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
const int kCacheVersion = 3;
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

int TPZKLRandomField::ChooseModes() const {
    const int n = (int)fLambda.size();
    if (fOpt.targetVarianceError > 0.) {
        for (int M = 1; M <= n; M++)
            if (VarianceError(M) <= fOpt.targetVarianceError) return M;
        return n;
    }
    return (fOpt.nModes > 0) ? std::min(fOpt.nModes, n) : n;
}

namespace {
// Formato do cache (versão 3; tudo por TPZBFileStream, binário nativo):
//   char[8] "SSR-KLC\n" | int versão, porder | double Lx, Ly, |Ω| | int64 neq, nel, nnós | uint64 assinatura da
//   malha KL | int64 nλ, M | double λ[nλ] | double Φ[neq x M] (coluna a coluna) | uint64 FNV-1a de tudo o que
//   precede (exceto o texto inicial) | char[8] "SSR-FIM\n"
// O tamanho do arquivo é conferido com o cabeçalho antes da leitura dos vetores (TPZBFileStream não acusa leitura
// além do fim) e o FNV-1a confere o conteúdo; o arquivo é gravado em <cache>.tmp<pid> e renomeado (atômico).
const char kCacheMagic[8] = {'S', 'S', 'R', '-', 'K', 'L', 'C', '\n'};
const char kCacheEnd[8] = {'S', 'S', 'R', '-', 'F', 'I', 'M', '\n'};
const uint64_t kCacheHeaderBytes = 8 + 2 * sizeof(int) + 3 * sizeof(double) + 3 * sizeof(int64_t) + sizeof(uint64_t) +
                                   2 * sizeof(int64_t);

struct TCacheHeader {
    int version = -1, porder = -1;
    double Lx = -1., Ly = -1., area = -1.;
    int64_t neq = -1, nel = -1, nnod = -1;
    uint64_t meshSig = 0;
    int64_t nl = -1, M = -1;
    uint64_t Hash() const {
        uint64_t h = fileio::Fnv1a("SSR-KLC", 7);
        const int64_t iv[8] = {version, porder, neq, nel, nnod, (int64_t)meshSig, nl, M};
        const double dv[3] = {Lx, Ly, area};
        h = fileio::HashWords(iv, 8, h);
        return fileio::HashWords(dv, 3, h);
    }
    uint64_t FileBytes() const {
        return kCacheHeaderBytes + sizeof(double) * ((uint64_t)nl + (uint64_t)neq * (uint64_t)M) + sizeof(uint64_t) + 8;
    }
};
} // namespace

bool TPZKLRandomField::ReadCache() {
    if (fOpt.cacheFile.empty()) return false;
    const std::string &file = fOpt.cacheFile;
    uint64_t fsize = 0;
    if (!fileio::FileSize(file, fsize)) return false;
    auto reject = [&](const std::string &why) {
        std::cout << "[KL] cache " << file << " " << why << "; recalculando\n";
        return false;
    };
    if (fsize < kCacheHeaderBytes) return reject("truncado (" + std::to_string(fsize) + " bytes)");
    TPZBFileStream in;
    in.OpenRead(file);
    if (!in.AmIOpenForRead()) return reject("ilegível");
    uint64_t magic = 0;
    TCacheHeader h;
    in.Read(&magic, 1);
    in.Read(&h.version, 1);
    in.Read(&h.porder, 1);
    in.Read(&h.Lx, 1);
    in.Read(&h.Ly, 1);
    in.Read(&h.area, 1);
    in.Read(&h.neq, 1);
    in.Read(&h.nel, 1);
    in.Read(&h.nnod, 1);
    in.Read(&h.meshSig, 1);
    in.Read(&h.nl, 1);
    in.Read(&h.M, 1);
    if (magic != fileio::Word8(kCacheMagic) || h.version != kCacheVersion)
        return reject("de outro formato ou versão");
    const int64_t neq = NEquations();
    const uint64_t sig = TSlopeGeometry::Signature(*fGMesh);
    if (h.porder != fOpt.porder || h.Lx != fOpt.Lx || h.Ly != fOpt.Ly || h.neq != neq ||
        h.nel != fGMesh->NElements() || h.nnod != fGMesh->NNodes() || h.meshSig != sig)
        return reject("incompatível (Lx, Ly, ordem ou malha KL diferentes)");
    if (h.nl < 1 || h.nl > neq || h.M < 1 || h.M > h.nl || !(h.area > 0.))
        return reject("corrompido (cabeçalho inválido)");
    if (fsize != h.FileBytes())
        return reject("com tamanho inesperado (" + std::to_string(fsize) + " bytes, esperado " +
                      std::to_string(h.FileBytes()) + "): truncado ou corrompido");
    std::vector<REAL> lambda(h.nl);
    in.Read(lambda.data(), (int)h.nl);
    TPZFMatrix<STATE> phi(neq, h.M);
    for (int64_t k = 0; k < h.M; k++) in.Read(&phi(0, k), (int)neq);
    uint64_t sum = 0, end = 0;
    in.Read(&sum, 1);
    in.Read(&end, 1);
    uint64_t check = fileio::HashWords(lambda.data(), lambda.size(), h.Hash());
    check = fileio::HashWords(&phi.g(0, 0), (size_t)(neq * h.M), check);
    if (sum != check || end != fileio::Word8(kCacheEnd)) return reject("corrompido (soma de verificação)");

    fLambda = lambda;
    fArea = h.area;
    const int Mreq = ChooseModes();
    if (Mreq > h.M) {
        fLambda.clear();
        return reject("com " + std::to_string(h.M) + " modos (pedidos " + std::to_string(Mreq) + ")");
    }
    fM = Mreq;
    if (fM == h.M) {
        fPhi = phi;
    } else {
        fPhi.Redim(neq, fM);
        for (int k = 0; k < fM; k++)
            for (int64_t i = 0; i < neq; i++) fPhi(i, k) = phi(i, k);
    }
    std::cout << "[KL] lido de " << file << ": " << fM << " de " << h.M << " modos gravados\n";
    return true;
}

void TPZKLRandomField::WriteCache() const {
    if (fOpt.cacheFile.empty()) return;
    TCacheHeader h;
    h.version = kCacheVersion;
    h.porder = fOpt.porder;
    h.Lx = fOpt.Lx;
    h.Ly = fOpt.Ly;
    h.area = fArea;
    h.neq = NEquations();
    h.nel = fGMesh->NElements();
    h.nnod = fGMesh->NNodes();
    h.meshSig = TSlopeGeometry::Signature(*fGMesh);
    h.nl = (int64_t)fLambda.size();
    h.M = fM;
    if (fPhi.Rows() != h.neq || fPhi.Cols() != h.M) return;
    uint64_t sum = fileio::HashWords(fLambda.data(), fLambda.size(), h.Hash());
    sum = fileio::HashWords(&fPhi.g(0, 0), (size_t)(h.neq * h.M), sum);
    const std::string tmp = fileio::TmpName(fOpt.cacheFile);
    {
        TPZBFileStream out;
        out.OpenWrite(tmp);
        if (!out.AmIOpenForWrite()) {
            std::cout << "[KL] não foi possível gravar o cache " << tmp << "\n";
            return;
        }
        const uint64_t magic = fileio::Word8(kCacheMagic), end = fileio::Word8(kCacheEnd);
        out.Write(&magic, 1);
        out.Write(&h.version, 1);
        out.Write(&h.porder, 1);
        out.Write(&h.Lx, 1);
        out.Write(&h.Ly, 1);
        out.Write(&h.area, 1);
        out.Write(&h.neq, 1);
        out.Write(&h.nel, 1);
        out.Write(&h.nnod, 1);
        out.Write(&h.meshSig, 1);
        out.Write(&h.nl, 1);
        out.Write(&h.M, 1);
        out.Write(fLambda.data(), (int)h.nl);
        for (int64_t k = 0; k < h.M; k++) out.Write(&fPhi.g(0, k), (int)h.neq);
        out.Write(&sum, 1);
        out.Write(&end, 1);
        out.CloseWrite();
    }
    uint64_t size = 0;
    if (!fileio::FileSize(tmp, size) || size != h.FileBytes()) {
        std::cout << "[KL] gravação incompleta do cache " << tmp << " (disco cheio?); cache não gravado\n";
        std::remove(tmp.c_str());
        return;
    }
    try {
        fileio::Commit(tmp, fOpt.cacheFile);
        std::cout << "[KL] cache gravado em " << fOpt.cacheFile << " (" << h.M << " modos)\n";
    } catch (std::exception &e) {
        std::cout << "[KL] cache não gravado: " << e.what() << "\n";
    }
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

    fM = ChooseModes();
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
