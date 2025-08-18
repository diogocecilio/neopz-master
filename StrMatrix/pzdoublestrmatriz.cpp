#include "pzdoublestrmatriz.h"
#include "TPZMatKLKernel.h"
#include "TPZStructMatrix.h"
#include "pzsubcmesh.h"
#include "TPZElementMatrixT.h"
#include "pzcompel.h"
#include "pzinterpolationspace.h"
#include "pzconnect.h"
#include "pzgeoel.h"
#include "tpzintpoints.h"
#include <cmath>

template<class TVar>
pzdoublestrmatriz<TVar>::pzdoublestrmatriz(TPZCompMesh *mesh)
: TPZStructMatrixT<TVar>(mesh)
{
}

template<class TVar>
pzdoublestrmatriz<TVar>::pzdoublestrmatriz(TPZAutoPointer<TPZCompMesh> cmesh)
: TPZStructMatrixT<TVar>(cmesh)
{
}

template<class TVar>
TPZMatrix<TVar>* pzdoublestrmatriz<TVar>::Create()
{
    return new TPZFMatrix<TVar>();
}

template<class TVar>
TPZBaseMatrix* pzdoublestrmatriz<TVar>::CreateAssemble(TPZBaseMatrix &rhs,
                                                       TPZAutoPointer<TPZGuiInterface> gui)
{
    // interface pede BaseMatrix; seguimos o fluxo padrão
    this->EquationFilter().Reset(); // protótipo: sem filtro
    rhs.Redim(0,0);

    this->InitCreateAssemble();

    auto *M = new TPZFMatrix<TVar>();
    // preenche M em Assemble(stiffness, rhs,...)
    Assemble(*M, rhs, gui);

    this->EndCreateAssemble(M);
    return M; // TPZEigenAnalysis fará dynamic_cast<TPZMatrix<TVar>*>(...)
}

template<class TVar>
void pzdoublestrmatriz<TVar>::Assemble(TPZBaseMatrix &stiffness, TPZBaseMatrix &rhs,
                                       TPZAutoPointer<TPZGuiInterface> /*gui*/)
{
    rhs.Redim(0,0);
    auto *M = dynamic_cast<TPZFMatrix<TVar>*>(&stiffness);
    if(!M) DebugStop();

    if (CurrentTargetIsA()) AssembleC(*M);
    else                    AssembleB(*M);
}

template<class TVar>
void pzdoublestrmatriz<TVar>::Assemble(TPZBaseMatrix &rhs,
                                       TPZAutoPointer<TPZGuiInterface> /*gui*/)
{
    rhs.Redim(0,0); // sem RHS neste problema
}

template<class TVar>
TPZStructMatrix* pzdoublestrmatriz<TVar>::Clone()
{
    return new pzdoublestrmatriz<TVar>(*this);
}

template<class TVar>
bool pzdoublestrmatriz<TVar>::CurrentTargetIsA() const
{
    // olha o primeiro material KL encontrado
    for (const auto &it : this->Mesh()->MaterialVec()) {
        if (auto *kl = dynamic_cast<TPZMatKLKernel*>(it.second)) {
            return kl->CurrentTarget() == TPZMatKLKernel::Target::A;
        }
    }
    return true; // fallback
}
template<class TVar>
void pzdoublestrmatriz<TVar>::AssembleC(TPZFMatrix<TVar> &C)
{
    TPZCompMesh* cmesh = this->Mesh();
    const int64_t neq = cmesh->NEquations();

    C.Redim(neq, neq);
    C.Zero();

    // --- 1) Filtra elementos válidos de VOLUME e guarda ponteiros úteis
    struct Entry {
        TPZInterpolationSpace* el = nullptr;
        TPZMatKLKernel*        mat = nullptr;
        int64_t                gelIndex = -1;
    };
    std::vector<Entry> vol; vol.reserve(cmesh->NElements());

    auto& elvec = cmesh->ElementVec();
    for (int64_t iel = 0; iel < (int64_t)cmesh->NElements(); ++iel) {
        auto *elx = dynamic_cast<TPZInterpolationSpace*>(elvec[iel]);
        if (!elx || !elx->Reference() || elx->HasDependency()) continue;

        // ignora BCs e materiais não-KL
        if (dynamic_cast<TPZBndCond*>(elx->Material())) continue;
        auto *matKL = dynamic_cast<TPZMatKLKernel*>(elx->Material());
        if (!matKL) continue;

        vol.push_back({elx, matKL, iel});
    }

    const std::size_t nvol = vol.size();
    if (nvol == 0) return;

    // --- 2) Pré-computa mapa local->global (dst) por elemento
    std::vector<TPZManVector<long>> dest; dest.resize(nvol);
    for (std::size_t a = 0; a < nvol; ++a) {
        // 'nx' deve bater com linhas do ce.fMat calculado para este elemento a
        const int nx = vol[a].el->NShapeF();
        TPZManVector<long> src, dstA;
        GetDestIndex(vol[a].gelIndex, nx, src, dstA);
        dest[a] = std::move(dstA);
    }

    // --- 3) Loop só no triângulo superior: b >= a (simetria)
    for (std::size_t a = 0; a < nvol; ++a) {
        for (std::size_t b = a; b < nvol; ++b) {

            TPZElementMatrixT<STATE> ce(cmesh, TPZElementMatrix::EK);
            try {
                vol[a].mat->CalcStiffNystrom(vol[a].el, vol[b].el, ce);
            } catch (...) {
                std::cout << "CalcStiffNystrom falhou: a="<<a<<" b="<<b
                << " (gel "<<vol[a].gelIndex<<","<<vol[b].gelIndex<<")\n";
                throw;
            }

            const int nx = ce.fMat.Rows();
            const int ny = ce.fMat.Cols();

            auto &IA = dest[a];
            auto &IB = dest[b];

            // segurança: pode acontecer de NShapeF() mudar após integração
            if ((int)IA.size() < nx) { TPZManVector<long> s,d; GetDestIndex(vol[a].gelIndex, nx, s, IA); }
            if ((int)IB.size() < ny) { TPZManVector<long> s,d; GetDestIndex(vol[b].gelIndex, ny, s, IB); }

            for (int i = 0; i < nx; ++i) {
                const long I = IA[i];
                for (int j = 0; j < ny; ++j) {
                    const long J = IB[j];
                    const TVar val = (TVar)ce.fMat(i,j);

                    // triângulo superior
                    C(I,J) += val;

                    if (b != a) {
                        // espelhamento (Hermitiano se TVar for complexo)
                        if constexpr (std::is_same_v<TVar,std::complex<float>> ||
                            std::is_same_v<TVar,std::complex<double>>) {
                            C(J,I) += std::conj(val);
                            } else {
                                C(J,I) += val;
                            }
                    }
                }
            }
        }
    }

    // --- 4) (Opcional) força simetria numérica
    // for (int64_t i = 0; i < neq; ++i) {
    //     for (int64_t j = i+1; j < neq; ++j) {
    //         TVar v = TVar(0.5)*(C(i,j) +
    //                 (std::is_complex_v<TVar> ? std::conj(C(j,i)) : C(j,i)));
    //         C(i,j) = v;
    //         C(j,i) = std::is_complex_v<TVar> ? std::conj(v) : v;
    //     }
    // }
}
/*
template<class TVar>
void pzdoublestrmatriz<TVar>::AssembleC(TPZFMatrix<TVar> &C)
{
    auto cmesh = this->Mesh();
    const long nelem = cmesh->NElements();
    auto &elvec = cmesh->ElementVec();
    const int64_t neq = cmesh->NEquations();

    C.Redim(neq, neq);
    C.Zero();

    const bool useGalerkin = (fCAssembly == ECAssembly::Galerkin);

    for (long iel=0; iel<nelem; ++iel) {
        auto *elx = dynamic_cast<TPZInterpolationSpace*>(elvec[iel]);
        if (!elx || !elx->Reference() || elx->HasDependency()) continue;

        auto *matKL = dynamic_cast<TPZMatKLKernel*>(elx->Material());
        if (!matKL) continue;

        for (long jel=0; jel<nelem; ++jel) {
            auto *ely = dynamic_cast<TPZInterpolationSpace*>(elvec[jel]);
            if (!ely || !ely->Reference() || ely->HasDependency()) continue;
            if (!dynamic_cast<TPZMatKLKernel*>(ely->Material()) ) continue;
            if(dynamic_cast<TPZBndCond*>(elx->Material())||dynamic_cast<TPZBndCond*>(ely->Material()))continue;
            if(dynamic_cast<TPZMatKLKernel*>(elx->Material())==nullptr||dynamic_cast<TPZMatKLKernel*>(ely->Material())==nullptr)continue;
            TPZElementMatrixT<STATE> ce(cmesh, TPZElementMatrix::EK);
            //elx->Material()->Print();
            auto *mat =
            dynamic_cast<TPZMatSingleSpace*>(elx->Material());
            if(!mat)
            {
                DebugStop();
            }
            try{
                if (useGalerkin) {
                    // NOVO: Galerkin padrão (duas regras de integração)
                    matKL->CalcStiffGalerkin(elx, ely, ce);
                } else {
                    // ANTIGO: Nyström
                    matKL->CalcStiffNystrom(elx, ely, ce);
                }
            } catch (...) {
                std::cout << "C-block falhou em iel=" << iel << " jel=" << jel << std::endl;
                throw;
            }

            const int nx = ce.fMat.Rows();
            const int ny = ce.fMat.Cols();

            TPZManVector<long> srcX, dstX, srcY, dstY;
            GetDestIndex(iel, nx, srcX, dstX);
            GetDestIndex(jel, ny, srcY, dstY);

            for (int i=0;i<nx;i++){
                const long I = dstX[i];
                for (int j=0;j<ny;j++){
                    const long J = dstY[j];
                    C(I,J) += ce.fMat(i,j);
                }
            }
        }
    }
}*/

template<class TVar>
void pzdoublestrmatriz<TVar>::AssembleB(TPZFMatrix<TVar> &B)
{
    auto cmesh = this->Mesh();
    const long nelem = cmesh->NElements();
    auto &elvec = cmesh->ElementVec();
    const int64_t neq = cmesh->NEquations();

    B.Redim(neq, neq);
    B.Zero();

    for (long iel=0; iel<nelem; ++iel) {
        auto *el = dynamic_cast<TPZInterpolationSpace*>(elvec[iel]);
        if (!el || !el->Reference() || el->HasDependency()) continue;

        auto *matKL = dynamic_cast<TPZMatKLKernel*>(el->Material());
        if (!matKL) continue;

        TPZElementMatrixT<STATE> be(cmesh, TPZElementMatrix::EK);
        matKL->CalcStiffMass(el, be, /*qmass=*/fQMass);

        const int n = be.fMat.Rows();
        TPZManVector<long> src, dst;
        GetDestIndex(iel, n, src, dst);

        for (int i=0;i<n;i++){
            const long I = dst[i];
            for (int j=0;j<n;j++){
                const long J = dst[j];
                B(I,J) += be.fMat(i,j);
            }
        }
    }
}

template<class TVar>
void pzdoublestrmatriz<TVar>::GetDestIndex(long iel, int nshape,
                                           TPZManVector<long> &source,
                                           TPZManVector<long> &dest)
{
    source.Resize(nshape);
    dest.Resize(nshape);

    long fullmatindex = 0L;
    long destindex = 0L;

    TPZCompEl *cel = this->Mesh()->ElementVec()[iel];
    const int ncon = cel->NConnects();
    for (int ic=0; ic<ncon; ++ic) {
        const long cindex = this->Mesh()->ElementVec()[iel]->ConnectIndex(ic);
        TPZConnect &c = this->Mesh()->ConnectVec()[cindex];

        const long seq = c.SequenceNumber();
        const long firsteq = this->Mesh()->Block().Position(seq);
        const int  ndf = this->Mesh()->Block().Size(seq);

        if (c.HasDependency() || c.IsCondensed()) {
            fullmatindex += ndf;
            continue;
        }
        for (int k=0; k<ndf; ++k) {
            source[destindex] = fullmatindex++;
            dest[destindex++] = firsteq + k;
        }
    }
}

template<class TVar>
int pzdoublestrmatriz<TVar>::ClassId() const
{
    // única para cada TVar
    return Hash("pzdoublestrmatriz") ^ (TPZStructMatrixT<TVar>::ClassId() << 1);
}

template<class TVar>
void pzdoublestrmatriz<TVar>::Write(TPZStream &buf, int withclassid) const
{
    TPZStructMatrixT<TVar>::Write(buf, withclassid);
    buf.Write(&fQMass, 1);
}

template<class TVar>
void pzdoublestrmatriz<TVar>::Read(TPZStream &buf, void *context)
{
    TPZStructMatrixT<TVar>::Read(buf, context);
    buf.Read(&fQMass, 1);
}

// ==== instanciacoes explicitas ====
template class pzdoublestrmatriz<STATE>;
template class pzdoublestrmatriz<CSTATE>;
