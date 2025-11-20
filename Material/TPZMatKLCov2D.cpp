#include "TPZMatKLCov2D.h"

#include "pzcmesh.h"
#include "pzintel.h"
#include "pzgeoel.h"
#include "tpzintpoints.h"
#include "TPZMaterialDataT.h"
#include "TPZBndCond.h"

#include <cmath>
#include <algorithm>

// ---------------- ctors ----------------
TPZMatKLCov2D::TPZMatKLCov2D(int id, int dim)
: TBase(), fDim(dim)
{
    this->SetId(id);
    // kernel default
    fKernel = [this](const TPZVec<REAL>& x, const TPZVec<REAL>& y)->STATE {
        return ExpKernel(x,y,fLx,fLy,fSigma2);
    };
}

// ------------- MASSA L² (B) -------------
void TPZMatKLCov2D::Contribute(const TPZMaterialDataT<STATE>& data,
                               REAL weight,
                               TPZFMatrix<STATE>& ek,
                               TPZFMatrix<STATE>& ef)
{
    const auto& phi = data.phi; // nshape x 1
    const int n = phi.Rows();
    const REAL detjac = data.detjac;

    for (int i = 0; i < n; i++) {
        const STATE phii = phi(i,0);
        for (int j = 0; j < n; j++) {
            ek(i,j) += phii * phi(j,0) * fabs(detjac) * weight;
        }
    }
    (void)ef;
}

// ------------- pós-processo -------------
int TPZMatKLCov2D::VariableIndex(const std::string &name) const
{
    if (name == "u" || name == "EVP_U") return 0;
    return -1;
}

int TPZMatKLCov2D::NSolutionVariables(int var) const
{
    if (var == 0) return 1;
    return 0;
}

void TPZMatKLCov2D::Solution(const TPZMaterialDataT<STATE>& data,
                             int var, TPZVec<STATE>& sol)
{
    if (var == 0) {
        sol.Resize(1);
        //sol[0] = (data.sol.empty() ? 0. : data.sol[0][0]);
    }
}
void TPZMatKLCov2D::Errors(const TPZMaterialDataT<STATE>& data,
                           TPZVec<double>& values)
{
    values.Resize(3);
    values.Fill(0.0);
    (void)data;
}
// ------------- serialização -------------
int TPZMatKLCov2D::ClassId() const
{
    return Hash("TPZMatKLCov2D") ^ (TBase::ClassId() << 1);
}

void TPZMatKLCov2D::Write(TPZStream& buf, int withclassid) const
{
    TBase::Write(buf, withclassid);
    buf.Write(&fDim,1);
    buf.Write(&fLx ,1);
    buf.Write(&fLy ,1);
    buf.Write(&fSigma2,1);
}

void TPZMatKLCov2D::Read(TPZStream& buf, void* context)
{
    TBase::Read(buf, context);
    buf.Read(&fDim,1);
    buf.Read(&fLx ,1);
    buf.Read(&fLy ,1);
    buf.Read(&fSigma2,1);
}

// ------------- kernel default -------------
STATE TPZMatKLCov2D::ExpKernel(const TPZVec<REAL>& x,
                               const TPZVec<REAL>& y,
                               REAL Lx, REAL Ly, STATE sigma2)
{
    const REAL dx = fabs(x[0]-y[0]);
    const REAL dy = fabs(x[1]-y[1]);
    return exp(-dx/Lx - dy/Ly);
}

// ------------- Nyström helpers -------------
void TPZMatKLCov2D::ApplyPhiColumn(int64_t j,
                                   const std::vector<GPEntry>& GP,
                                   std::vector<double>& t)
{
    t.assign(GP.size(), 0.0);
    for (size_t m=0; m<GP.size(); ++m) {
        const auto &g = GP[m];
        for (size_t k=0; k<g.dof.size(); ++k) {
            if (g.dof[k] == j) { t[m] += g.phi[k]; }
        }
    }
}

void TPZMatKLCov2D::PhiT_times_vec(const std::vector<GPEntry>& GP,
                                   const std::vector<double>& s,
                                   TPZFMatrix<STATE>& y)
{
    for (size_t m=0; m<GP.size(); ++m) {
        const auto &g = GP[m]; const double sm = s[m];
        for (size_t k=0; k<g.dof.size(); ++k) {
            const int64_t i = g.dof[k];
            y(i,0) += g.phi[k] * sm;
        }
    }
}

void TPZMatKLCov2D::KernelApply(const std::vector<GPEntry>& GP,
                                const KernelType& ker,
                                const std::vector<double>& r,
                                std::vector<double>& s)
{
    const size_t M = GP.size();
    s.assign(M, 0.0);
    for (size_t m=0; m<M; ++m) {
        double acc = 0.0;
        const auto &xm = GP[m].x;
        for (size_t n=0; n<M; ++n) {
            acc += ker(xm, GP[n].x) * r[n];
        }
        s[m] = acc;
    }
}

#include <cmath>
#include <algorithm>
#include "pzconnect.h" // <- necessário para con.Order()

void TPZMatKLCov2D::BuildGlobalGPCatalog(TPZCompMesh& cmesh,
                                         int qorder,
                                         std::vector<GPEntry>& GP,
                                         int64_t& neq)
{
    GP.clear();
    neq = cmesh.NEquations();

    TPZGeoMesh* gmesh = cmesh.Reference();
    const int dim = gmesh ? gmesh->Dimension() : cmesh.Dimension();

    const int64_t nel = cmesh.NElements();
    for (int64_t icel = 0; icel < nel; ++icel)
    {
        TPZCompEl* cel = cmesh.Element(icel);
        auto* intel = dynamic_cast<TPZInterpolationSpace*>(cel);
        if (!intel) continue;

        TPZGeoEl* gel = intel->Reference();
        if (!gel || gel->Dimension() != dim) continue;

        // regra de integração do lado interior
        const int side = gel->NSides() - 1;
        TPZAutoPointer<TPZIntPoints> ir = gel->CreateSideIntegrationRule(side, qorder);
        const int npts = ir->NPoints();

        // (opcional) data, caso queira evoluir para usar ComputeRequiredData
        TPZMaterialDataT<STATE> data;
        intel->InitMaterialData(data);

        for (int ip = 0; ip < npts; ++ip)
        {
            // --- ponto e peso no ref ---
            TPZManVector<REAL,3> qsi(3,0.0), x(3,0.0);
            REAL wloc = 0.0;
            ir->Point(ip, qsi, wloc);

            // --- mapeia para físico e |detJ| ---
            TPZFMatrix<REAL> J, axes, jinv;

            intel->ComputeRequiredData(data, qsi);

            // coords físicas e peso com |detJ|
            GPEntry gp;
            gp.x = { data.x[0], data.x[1] };
            gp.w = wloc * std::abs(data.detjac);

            gp.dof.clear();
            gp.phi.clear();

            // varre connects e associa as φ aos DOFs na MESMA ordem
            const int ncon = intel->NConnects();
            int offset = 0; // deslocamento na lista de shapes
            for (int ic = 0; ic < ncon; ++ic)
            {
                // pega o connect e a ordem a partir dele
                const int cindex = intel->ConnectIndex(ic);
                if (cindex < 0) {
                    // Connect inválido: não deve acontecer em H1 com malha bem construída;
                    // por segurança, não associamos DOFs e tentamos avançar usando ordem 1.
                    const int nshc_fallback = intel->NConnectShapeF(ic, /*ord=*/1);
                    offset += nshc_fallback;
                    continue;
                }

                TPZConnect& con = cmesh.ConnectVec()[cindex];
                const int ord = con.Order();                       // <- substitui ConnectOrder(ic)
                const int nshc = intel->NConnectShapeF(ic, ord);   // #shapes deste connect

                const int seq = con.SequenceNumber();
                if (seq >= 0) {
                    const int64_t eqstart = cmesh.Block().Position(seq);
                    const int blksz = cmesh.Block().Size(seq); // #DOFs nesse connect

                    // Em H1 escalar, blksz == nshc. Se não for, emparelha pela menor.
                    const int take = std::min(blksz, nshc);
                    for (int s = 0; s < take; ++s) {
                        gp.dof.push_back(eqstart + s);
                        gp.phi.push_back( data.phi(offset + s, 0) );
                    }
                }

                offset += nshc; // avança para o próximo bloco de shapes
            }

            // sanidade
            if (gp.dof.size() != gp.phi.size()) {
                std::cout << "[WARN] gp.dof.size()=" << gp.dof.size()
                << " != gp.phi.size()=" << gp.phi.size()
                << "  (icel=" << icel << ", ip=" << ip << ")\n";
            }

            GP.emplace_back(std::move(gp));
        }
    }

    // mantém neq coerente
    neq = cmesh.NEquations();
}

// ------------- Nyström explícito -------------
void TPZMatKLCov2D::BuildB_Nystrom(TPZCompMesh& cmesh, int qorder,
                                   TPZFMatrix<STATE>& B)
{
    std::vector<GPEntry> GP; int64_t neq = 0;
    BuildGlobalGPCatalog(cmesh, qorder, GP, neq);

    // use a dimensão oficial da malha
    const int64_t N = neq;

    // (opcional) sanidade: verifique se todos ids estão no intervalo [0, N-1]
    #ifdef DEBUG
    for (size_t m = 0; m < GP.size(); ++m) {
        for (auto id : GP[m].dof) {
            if (id < 0 || id >= N) {
                std::cout << "[WARN] GP " << m << " id=" << id
                << " fora de [0," << (N-1) << "]\n";
            }
        }
    }
    #endif

    B.Redim(N, N);
    B.Zero();

    std::vector<double> t(GP.size());
    for (int64_t j=0; j<N; ++j) {
        ApplyPhiColumn(j, GP, t);                           // t = Φ e_j
        for (size_t m=0; m<GP.size(); ++m) t[m] *= GP[m].w; // W t
        TPZFMatrix<STATE> y(N,1); y.Zero();
        PhiT_times_vec(GP, t, y);                           // y = Φᵀ (W t)
        for (int64_t i=0;i<N;i++) B(i,j) = y(i,0);
    }
    for (int64_t i=0;i<N;i++)
        for (int64_t j=i+1;j<N;j++){
            const double a = 0.5*(B(i,j)+B(j,i));
            B(i,j)=B(j,i)=a;
        }
}

void TPZMatKLCov2D::BuildC_Nystrom(TPZCompMesh& cmesh, int qorder,
                                   const KernelType& ker,
                                   TPZFMatrix<STATE>& C)
{
    // ====== CONTROLES DE IMPRESSÃO ======
    const int  TRACE_EQJ      = 3;   // coluna (eq) que você quer rastrear; -1 = não rastrear
    const int  TRACE_M_LIMIT  = 16;  // quantos GPs imprimir (t/s) nos vetores
    const bool TRACE_GPS      = false;   // imprimir resumo dos GPs (x,w, ndof) no início
    const bool TRACE_VECS     = false;   // imprimir t, Wt, s, Ws para a coluna rastreada
    const bool TRACE_COL_Y    = false;   // imprimir a coluna final y (= C(:,j)) na ordem de eq

    auto print_vec = [&](const char* name, const std::vector<double>& v){
        std::cout << name << " (primeiros " << std::min<int>(TRACE_M_LIMIT,(int)v.size())
        << " de " << v.size() << ")\n";
        std::cout << std::scientific << std::setprecision(16);
        for (int i = 0; i < (int)v.size() && i < TRACE_M_LIMIT; ++i){
            std::cout << i << "  " << v[i] << "\n";
        }
    };

    auto print_col_y = [&](const char* name, const TPZFMatrix<STATE>& y, int j){
        std::cout << name << "  C(:, " << j << ")\n";
        std::cout << std::scientific << std::setprecision(16);
        for (int i = 0; i < y.Rows(); ++i){
            std::cout << i << "  " << y.GetVal(i,0) << "\n";
        }
    };

    // ====== CATÁLOGO GLOBAL DE GPS ======
    std::vector<GPEntry> GP; int64_t neq = 0;
    BuildGlobalGPCatalog(cmesh, qorder, GP, neq);

    int64_t N = 0;
    for (auto& g: GP) for (auto id: g.dof) N = std::max<int64_t>(N, id+1);
    C.Redim(N,N); C.Zero();

    // — Prints de sanidade dos GPs —
    if (TRACE_GPS){
        double sumw = 0.0;
        for (auto& g : GP) sumw += g.w;
        std::cout << "[KL] GP.size = " << GP.size()
        << " | N(eq) = " << N
        << " | sum(w) = " << std::scientific << std::setprecision(16) << sumw << "\n";

        int mlim = std::min<int>(TRACE_M_LIMIT, (int)GP.size());
        for (int m = 0; m < mlim; ++m){
            // Ajuste 'g.x' -> 'g.co' se seu GPEntry usa outro nome
            const auto& g = GP[m];
            double x0 = 0.0, x1 = 0.0;
            if (g.x.size() >= 2) { x0 = g.x[0]; x1 = g.x[1]; } // coords físicas
            std::cout << "GP " << m
            << "  x=(" << std::scientific << std::setprecision(16)
            << x0 << "," << x1 << ")"
            << "  w=" << g.w
            << "  ndof=" << g.dof.size() << "\n";
        }

        // Alguns k(xa,xb) de amostra
        if (GP.size() >= 2){
            std::cout << "k(GP0,GP0)=" << ker(GP[0].x, GP[0].x) << "\n";
            std::cout << "k(GP0,GP1)=" << ker(GP[0].x, GP[1].x) << "\n";
            int a = std::min<int>( (int)GP.size()-1, 10 );
            std::cout << "k(GP" << a << ",GP0)=" << ker(GP[a].x, GP[0].x) << "\n";
        }
    }

    // ====== LOOP DAS COLUNAS ======
    std::vector<double> t, s;
    for (int64_t j = 0; j < N; ++j) {

        // t = Φ e_j
        ApplyPhiColumn(j, GP, t);
        if (TRACE_VECS && (TRACE_EQJ == (int)j)) {
            print_vec("t = Phi * e_j", t);
        }

        // W * t
        for (size_t m = 0; m < GP.size(); ++m) t[m] *= GP[m].w;
        if (TRACE_VECS && (TRACE_EQJ == (int)j)) {
            print_vec("Wt = W * t", t);
        }

        // s = K (W t)
        KernelApply(GP, ker, t, s);
        if (TRACE_VECS && (TRACE_EQJ == (int)j)) {
            print_vec("s = K * (W t)", s);
        }

        // W * s
        for (size_t m = 0; m < GP.size(); ++m) s[m] *= GP[m].w;
        if (TRACE_VECS && (TRACE_EQJ == (int)j)) {
            print_vec("Ws = W * s", s);
        }

        // y = Φᵀ (W s)  -> coluna j de C
        TPZFMatrix<STATE> y(N,1); y.Zero();
        PhiT_times_vec(GP, s, y);
        if (TRACE_COL_Y && (TRACE_EQJ == (int)j)) {
            print_col_y("y = Phi^T * (W s)  =>", y, (int)j);
        }

        for (int64_t i = 0; i < N; i++) C(i,j) = y(i,0);
    }

    // simetrização (como no seu código)
    for (int64_t i = 0; i < N; i++)
        for (int64_t j = i+1; j < N; j++) {
            const double a = 0.5*(C(i,j) + C(j,i));
            C(i,j) = C(j,i) = a;
        }
}
