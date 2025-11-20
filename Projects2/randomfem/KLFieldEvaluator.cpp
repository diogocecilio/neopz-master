
// ================================
// File: KLFieldEvaluator.cpp  (versão em classe)
// ================================
#include "KLFieldEvaluator.h"
#include "pzconnect.h"
#include "pzcompel.h"

KLFieldEvaluator::KLFieldEvaluator(const TPZCompMesh* base_cmesh,
                                   const TPZFMatrix<STATE>& PHI_nodal,
                                   const TPZVec<STATE>& lambdas)
: fBaseCmesh(base_cmesh), fPHI(PHI_nodal), fL(lambdas), fM(PHI_nodal.Cols())
{
    if (!fBaseCmesh) DebugStop();
    if (fM <= 0 || fPHI.Rows() <= 0 || fL.size() != fM) DebugStop();
}

TPZFMatrix<STATE> KLFieldEvaluator::DrawAlpha(int samples, uint64_t seed) const {
    TPZFMatrix<STATE> Alpha(fM, samples, 0.0);
    std::mt19937_64 rng(seed);
    std::normal_distribution<double> N01(0.0, 1.0);
    for (int k=0;k<fM;k++) {
        const STATE s = std::sqrt(fL[k]);
        for (int sidx=0;sidx<samples;sidx++) Alpha(k,sidx) = s * (STATE)N01(rng);
    }
    return Alpha;
}

void KLFieldEvaluator::SaveAlphaBIN(const std::string& path,
                                    const TPZFMatrix<STATE>& Alpha,
                                    const KLBaseManifest& man)
{
    std::ofstream f(path, std::ios::binary);
    if(!f) throw std::runtime_error("SaveAlphaBIN: não abriu arquivo");
    f.write((const char*)&man, sizeof(KLBaseManifest));
    const int M = Alpha.Rows(); const int S = Alpha.Cols();
    f.write((const char*)&M, sizeof(int));
    f.write((const char*)&S, sizeof(int));
    for (int i=0;i<M;i++) for (int j=0;j<S;j++) {
        STATE v = Alpha.GetVal(i,j); f.write((const char*)&v, sizeof(STATE));
    }
}

bool KLFieldEvaluator::LoadAlphaBIN(const std::string& path,
                                    TPZFMatrix<STATE>& Alpha,
                                    KLBaseManifest& man)
{
    std::ifstream f(path, std::ios::binary);
    if(!f) return false;
    f.read((char*)&man, sizeof(KLBaseManifest));
    int M=0,S=0; f.read((char*)&M,sizeof(int)); f.read((char*)&S,sizeof(int));
    Alpha.Redim(M,S);
    for (int i=0;i<M;i++) for (int j=0;j<S;j++) {
        STATE v; f.read((char*)&v,sizeof(STATE)); Alpha(i,j)=v;
    }
    return true;
}

// KLFieldEvaluator.cpp

#include "pzconnect.h"  // (garanta esse include)

void KLFieldEvaluator::ComputeShapeGlobalDOFIndices(
    TPZInterpolationSpace* intel,
    TPZCompMesh const* cmesh,
    TPZVec<long>& dofIdx)
{
    const int nshape = intel->NShapeF();
    dofIdx.Resize(nshape);

    int ish = 0;
    const int ncon = intel->NConnects();
    for (int ic=0; ic<ncon; ic++) {
        TPZConnect& C = intel->Connect(ic);
        const int seq = C.SequenceNumber();
        const int pos = cmesh->Block().Position(seq);

        int order = 0;
        // Use o que existir na sua branch:
        order = C.Order();                 // << preferido
        // se não compilar, troque por:
        // order = intel->ConnectOrder(ic);

        const int nshape_c = intel->NConnectShapeF(ic, order);
        for (int k=0; k<nshape_c; k++) {
            dofIdx[ish++] = pos + k;
        }
    }
}


// --- antes: static inline void EvalPhiK_at_qsi_on_base(...)
// --- depois:
void KLFieldEvaluator::EvalPhiK_at_qsi_on_base(
    TPZInterpolationSpace* intel_base,
    TPZCompMesh const* /*base_cmesh*/,
    const TPZFMatrix<STATE>& PHI_nodal,
    const TPZManVector<REAL,3>& qsi,
    TPZVec<long> const& dofIdx,
    TPZVec<STATE>& phiK_out)
{
    const int nshape = intel_base->NShapeF();
    const int M = PHI_nodal.Cols();
    phiK_out.Resize(M);
    for (int k=0;k<M;k++) phiK_out[k]=0.0;

    TPZFNMatrix<40,REAL>  phi;
    TPZFNMatrix<120,REAL> dphi;

    TPZManVector<REAL,3> qsi_loc = qsi;  // Shape pede qsi mutável
    intel_base->Shape(qsi_loc, phi, dphi);

    for (int a=0; a<nshape; a++) {
        const long g = dofIdx[a];
        const REAL Na = phi(a,0);
        for (int k=0;k<M;k++) {
            phiK_out[k] += Na * PHI_nodal.GetVal(g,k);
        }
    }
}


void KLFieldEvaluator::BuildPhiTableAtTargetQPs(TPZCompMesh* target_cmesh,
                                                TPZMatWithMem<TKLPointMem>* target_mat,
                                                TPZFMatrix<float>& PhiQP_out,
                                                TPZVec<STATE>& varG_out) const
{
    auto base_g = fBaseCmesh->Reference();
    auto tgt_g  = target_cmesh->Reference();
    base_g->BuildConnectivity();
    tgt_g->BuildConnectivity();
    target_cmesh->LoadReferences();

    long maxMem = -1;
    for (auto cel : target_cmesh->ElementVec()) {
        auto intel = dynamic_cast<TPZInterpolationSpace*>(cel);
        if (!intel) continue;
        if (!intel->Material() ||
            intel->Material() != dynamic_cast<TPZMaterial*>(target_mat)) continue;// <---
        TPZMaterialDataT<STATE> data; // <---
        intel->InitMaterialData(data);
        const int nq = intel->GetIntegrationRule().NPoints();
        for (int iq=0;iq<nq;iq++) {
            REAL w; TPZManVector<REAL,3> qsi(3,0.);
            intel->GetIntegrationRule().Point(iq, qsi, w);
            intel->ComputeRequiredData(data, qsi); // <---
            if (data.intGlobPtIndex > maxMem) maxMem = data.intGlobPtIndex;
        }
    }
    const long Nmem = maxMem+1; const int M = fM;
    PhiQP_out.Redim(Nmem, M);
    varG_out.Resize(Nmem); for (long i=0;i<Nmem;i++) varG_out[i]=0.;

    TPZVec<long> dofIdx; TPZVec<STATE> phiK; TPZManVector<REAL,3> qsi_base(3,0.);

    for (auto cel : target_cmesh->ElementVec()) {
        auto intel_tgt = dynamic_cast<TPZInterpolationSpace*>(cel);
        if (!intel_tgt) continue;
        if (!intel_tgt->Material() || intel_tgt->Material() != dynamic_cast<TPZMaterial*>(target_mat)) continue;
        TPZMaterialDataT<STATE> data; // <---
        intel_tgt->InitMaterialData(data);
        const auto& intrule = intel_tgt->GetIntegrationRule();
        const int nq = intrule.NPoints();
        for (int iq=0;iq<nq;iq++) {
            REAL w; TPZManVector<REAL,3> qsi_tgt(3,0.);
            intrule.Point(iq, qsi_tgt, w);
            data.intLocPtIndex = iq;
            intel_tgt->ComputeRequiredData(data, qsi_tgt); // <---
            const long memIdx = data.intGlobPtIndex;
            TPZManVector<REAL,3> xphys = data.x;

            long elid=-1; TPZGeoEl* gel_base = base_g->FindElement(xphys, qsi_base, elid, 2);
            if (!gel_base) continue;
            auto cel_base = dynamic_cast<TPZInterpolationSpace*>(gel_base->Reference());
            if (!cel_base) continue;

            ComputeShapeGlobalDOFIndices(cel_base, fBaseCmesh, dofIdx);
            EvalPhiK_at_qsi_on_base(cel_base, fBaseCmesh, fPHI, qsi_base, dofIdx, phiK);

            for (int k=0;k<M;k++) {
                const float v = (float)phiK[k];
                PhiQP_out(memIdx,k) = v;
                varG_out[memIdx] += fL[k] * (STATE)v * (STATE)v;
            }
        }
    }
}


void KLFieldEvaluator::BuildQPMap(TPZCompMesh* target_cmesh,
                                  TPZMatWithMem<TKLPointMem>* target_mat,
                                  std::vector<KLQPMapEntry>& out_map) const
{
    auto tgt_g = target_cmesh->Reference();
    auto base_g= fBaseCmesh->Reference();
    tgt_g->BuildConnectivity();
    base_g->BuildConnectivity();
    target_cmesh->LoadReferences();

    long maxMem=-1;
    for (auto cel : target_cmesh->ElementVec()) {
        auto intel = dynamic_cast<TPZInterpolationSpace*>(cel);
        if (!intel) continue;
        if (!intel->Material() || intel->Material() != dynamic_cast<TPZMaterial*>(target_mat)) continue;
        TPZMaterialDataT<STATE> data; // <---
        intel->InitMaterialData(data);
        const int nq = intel->GetIntegrationRule().NPoints();
        for (int iq=0;iq<nq;iq++) {
            REAL w; TPZManVector<REAL,3> qsi(3,0.);
            intel->GetIntegrationRule().Point(iq, qsi, w);
            intel->ComputeRequiredData(data, qsi); // <---
            if (data.intGlobPtIndex > maxMem) maxMem = data.intGlobPtIndex;
        }
    }
    const long Nmem = maxMem+1;
    out_map.clear(); out_map.resize(Nmem);

    for (auto cel : target_cmesh->ElementVec()) {
        auto intel = dynamic_cast<TPZInterpolationSpace*>(cel);
        if (!intel) continue;
        if (!intel->Material() || intel->Material() != dynamic_cast<TPZMaterial*>(target_mat)) continue;
        TPZMaterialDataT<STATE> data; // <---
        intel->InitMaterialData(data);
        const auto& intrule = intel->GetIntegrationRule();
        const int nq = intrule.NPoints();
        for (int iq=0;iq<nq;iq++) {
            REAL w; TPZManVector<REAL,3> qsi_t(3,0.);
            intrule.Point(iq, qsi_t, w);
            data.intLocPtIndex = iq;
            intel->ComputeRequiredData(data, qsi_t); // <---
            long memIdx = data.intGlobPtIndex;

            TPZManVector<REAL,3> qsi_b(3,0.), x = data.x;
            long elid=-1; TPZGeoEl* gel_base = fBaseCmesh->Reference()->FindElement(x, qsi_b, elid, 2);
            if (!gel_base) continue;
            out_map[memIdx].gel_id = gel_base->Index();
            out_map[memIdx].qsi = qsi_b;
        }
    }
}


void KLFieldEvaluator::BuildPhiTableFromQPMap(const std::vector<KLQPMapEntry>& qpmap,
                                              TPZFMatrix<float>& PhiQP_out,
                                              TPZVec<STATE>& varG_out) const
{
    const long Nmem = (long)qpmap.size();
    const int M = fM;
    PhiQP_out.Redim(Nmem, M);
    varG_out.Resize(Nmem); for (long i=0;i<Nmem;i++) varG_out[i]=0.;

    auto base_g = fBaseCmesh->Reference();
    TPZVec<long> dofIdx; TPZVec<STATE> phiK;

    for (long i=0;i<Nmem;i++) {
        const auto& e = qpmap[i]; if (e.gel_id < 0) continue;
        TPZGeoEl* gel_base = base_g->Element(e.gel_id); if (!gel_base) continue;
        auto cel_base = dynamic_cast<TPZInterpolationSpace*>(gel_base->Reference());
        if (!cel_base) continue;
        ComputeShapeGlobalDOFIndices(cel_base, fBaseCmesh, dofIdx);
        EvalPhiK_at_qsi_on_base(cel_base, fBaseCmesh, fPHI, e.qsi, dofIdx, phiK);
        for (int k=0;k<M;k++) {
            const float v = (float)phiK[k];
            PhiQP_out(i,k) = v;
            varG_out[i] += fL[k] * (STATE)v * (STATE)v;
        }
    }
}

void KLFieldEvaluator::ApplySampleToMemory(const TPZFMatrix<float>& PhiQP,
                                           const TPZVec<STATE>& varG,
                                           const TPZFMatrix<STATE>& Alpha,
                                           int sidx,
                                           bool use_lognormal,
                                           STATE mean_mu,
                                           STATE cov_c,
                                           TPZMatWithMem<TKLPointMem>* target_mat) const
{
    auto mem = target_mat->GetMemory(); // shared_ptr/TPZAutoPointer
    target_mat->SetUpdateMem(true);     // garantir escrita

    const long Nmem = PhiQP.Rows();
    const int  M    = PhiQP.Cols();

    STATE xi=0., lam=0.;
    if (use_lognormal) {
        const STATE c = cov_c; // σ/μ
        xi  = std::sqrt(std::log(1.0 + c*c));
        lam = std::log(mean_mu) - 0.5*xi*xi;
    }

    for (long i=0;i<Nmem;i++) {
        double g = 0.0;
        for (int k=0;k<M;k++) g += (double)PhiQP(i,k) * (double)Alpha.GetVal(k,sidx);
        double E = 0.0;
        if (use_lognormal) {
            const double inv = (varG[i] > 0) ? 1.0/std::sqrt((double)varG[i]) : 0.0;
            const double ghat = g * inv;
            E = std::exp((double)lam + (double)xi * ghat);
        } else {
            E = (double)mean_mu + (double)(cov_c*mean_mu) * g; // gaussiano
        }
        (*mem)[i].fE = (STATE)E; // <--- aqui
    }

    target_mat->SetUpdateMem(false);
}
