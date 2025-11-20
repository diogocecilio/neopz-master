// main.cpp
#include <random>
#include <cmath>
#include <string>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <set>
#include "Plasticity/TPZElasticResponse.h"
#include "pzgmesh.h"
#include "pzcmesh.h"
#include "pzgeoquad.h"
#include "TPZGeoLinear.h"
#include "tpzgeoelrefpattern.h"
#include "pzintel.h"
#include "pzstack.h"

#include "TPZEigenAnalysis.h"
#include "TPZKrylovEigenSolver.h"
#include "pzdoublestrmatriz.h"
#include "pzskylstrmatrix.h"

#include "TPZMatKLKernel.h"
#include "Elasticity/TPZMatElastic2DMem.h"
#include "Elasticity/TPZElasticMem.h"

#include "pzpostprocanalysis.h"
#include "pzstepsolver.h"
#include "TPZBFileStream.h"
#include "TPZVTKGeoMesh.h"
#include "Projection/TPZL2Projection.h"

// ------------------------------------------------------------
// Assinaturas (mantendo nomes)
// ------------------------------------------------------------
TPZGeoMesh * CreateGeoMeshMathematicaLike(int matId, int nRef);
TPZCompMesh * CompMeshElastic(TPZGeoMesh * gmesh, int porder, int matId);
TPZCompMesh * BuildCompMeshKL(TPZGeoMesh *gmesh, int porder, int matId);
TPZFMatrix<STATE> BuildPhiSqrtLambda(TPZCompMesh* cmesh,const TPZFMatrix<CSTATE>& eigenvectors,const TPZVec<CSTATE>&eigenvalues,int M);
void PostElastic(TPZCompMesh* cmesh,const std::string& vtkfile,int matid = 1);
void CreatePostProcessingMesh(TPZCompMesh* cmesh,TPZPostProcAnalysis* pproc,int matid);
void PostProcessVariables(TPZStack<std::string>& scal, TPZStack<std::string>& vec);
REAL UyAtPoint(TPZCompMesh* cmesh, REAL x, REAL y);
void ComputeField(TPZCompMesh* source,TPZCompMesh* target,TPZFMatrix<REAL> THETA,TPZFMatrix<CSTATE> evecs_rows,TPZVec<CSTATE> evals);

// ------------------------------------------------------------
// Utilitário: lê u_y em um ponto físico (via postproc de elemento)
// ------------------------------------------------------------
REAL UyAtPoint(TPZCompMesh* cmesh, REAL x, REAL y)
{
    TPZGeoMesh* gmesh = cmesh->Reference();
    TPZManVector<REAL,3> X(3,0.); X[0]=x; X[1]=y;
    TPZManVector<REAL,3> qsi(3,0.);
    int64_t elid=0;

    cmesh->LoadReferences();          // liga cada TPZGeoEl ao seu TPZCompEl
    TPZGeoEl* gel = gmesh->FindElement(X, qsi, elid, /*dim=*/2);
    if(!gel || !gel->Reference()) DebugStop();

    auto* cel = dynamic_cast<TPZInterpolationSpace*>(gel->Reference());
    if(!cel) DebugStop();

    int var=1;
    TPZVec<REAL> sol;
    cel->Solution(qsi, var, sol);
    return sol[1]; // componente y
}


TPZFMatrix<REAL> ComputeTheta ( int M, REAL mean, REAL cov,int samples)
{
    std::normal_distribution<double> distribution(0.,1.);
    std::mt19937 generator(std::random_device{}());

    TPZFMatrix<REAL>  THETA ( M, samples, 0. );
    for ( int isample = 0; isample < samples; isample++ )
        for ( int irdvar = 0; irdvar < M; irdvar++ )
            THETA ( irdvar,isample ) = distribution ( generator );

    return THETA;
}
// ------------------------------------------------------------
// PHI_sqrtLambda (ndof × nm). Mantido por compatibilidade.
// ------------------------------------------------------------
TPZFMatrix<STATE> BuildPhiSqrtLambda(
    TPZCompMesh* cmesh,
    const TPZFMatrix<CSTATE>& eigenvectors, // cols = modos
    const TPZVec<CSTATE>&     eigenvalues,  // tamanho >= nm
    int M)
{
      TPZFMatrix<STATE>soltemp=cmesh->Solution();
      const int ndof = cmesh->NEquations();
      const int nm   = std::min(M, (int)eigenvectors.Cols());
      //if (nm <= 0) return;

      TPZFMatrix<STATE> PHI_sqrtLambda(ndof, nm, (STATE)0);

      // materiais de volume (sem BC) para integração
      std::set<int> mats;
      for (auto &it : cmesh->MaterialVec()) {
        if (!it.second) continue;
        if (dynamic_cast<TPZBndCond*>(it.second)) continue;
        mats.insert(it.first);
      }

      TPZFMatrix<STATE> sol(ndof,1,(STATE)0);

      for (int k = 0; k < nm; ++k) {
        // 1) autovetor k (parte real)
        for (int i=0; i<ndof; i++) {
          sol(i,0) = (STATE)eigenvectors.GetVal(i,k).real();
        }

        // 2) normalização na métrica de massa: ||phi||_M = sqrt( ∫ phi^2 dΩ )
        double scale = 1.0;

        cmesh->LoadSolution(sol);

        TPZVec<STATE> Ivec = cmesh->Integrate("SolutionSquared", mats);
        const double I = (Ivec.size() ? (double)Ivec[0] : 0.0);
        if (I > 1e-30) scale = 1.0/std::sqrt(I);

        // 3) √λ_k (parte real, truncada em 0)
        const double lamk    = std::max(0.0, (double)eigenvalues[k].real());
        const double rootlam = std::sqrt(lamk);

        // 4) coluna k := √λ_k * (phi_k * scale)
        for (int i=0; i<ndof; i++) {
          PHI_sqrtLambda(i,k) = (STATE)(rootlam* (double)sol(i,0) * scale);
        }
      }

      cmesh->Solution().Zero();
      return PHI_sqrtLambda;
}

// ------------------------------------------------------------
// Pós-processo
// ------------------------------------------------------------
void PostProcessVariables(TPZStack<std::string>& scal, TPZStack<std::string>& vec)
{
    scal.Push("Young");
    scal.Push("Poisson");
    scal.Push("POrder");
    vec.Push("Displacement");
    vec.Push("Strain");
    vec.Push("Stress");
}

void CreatePostProcessingMesh(TPZCompMesh* cmesh,TPZPostProcAnalysis* pproc,int matid)
{
    if (pproc->ReferenceCompMesh() != cmesh) {
        pproc->SetCompMesh(cmesh);
        TPZStack<std::string> scal, vec, all;
        PostProcessVariables(scal, vec);
        for (auto i=0; i<scal.size();  ++i) all.Push(scal[i]);
        for (auto i=0; i<vec.size();   ++i) all.Push(vec[i]);
        TPZVec<int> matids(1); matids[0] = matid;
        pproc->SetPostProcessVariables(matids, all);
        TPZFStructMatrix<REAL> str(pproc->Mesh());
        str.SetNumThreads(0);
        pproc->SetStructuralMatrix(str);
    }
    pproc->TransferSolution();
}

void PostElastic(TPZCompMesh* cmesh,const std::string& vtkfile,int matid)
{
    TPZPostProcAnalysis pproc;
    CreatePostProcessingMesh(cmesh, &pproc, matid);
    TPZStack<std::string> scal, vec;
    PostProcessVariables(scal, vec);
    pproc.DefineGraphMesh(/*dim=*/2, scal, vec, vtkfile);
    pproc.PostProcess(0);
}

// ------------------------------------------------------------
// Geometria: quadrado [-0.5,0.5]^2 com 2 lados nomeados para BCs
// ------------------------------------------------------------
TPZGeoMesh * CreateGeoMeshMathematicaLike(int matId, int nRef)
{
    auto * gmesh = new TPZGeoMesh;
    gmesh->SetDimension(2);

    const double X[4][2] = { {-0.5,-0.5},{0.5,-0.5},{0.5,0.5},{-0.5,0.5} };
    gmesh->NodeVec().Resize(4);
    for (int i=0;i<4;i++){
        TPZManVector<REAL,3> c(3,0.); c[0]=X[i][0]; c[1]=X[i][1];
        gmesh->NodeVec()[i].Initialize(c, *gmesh);
    }

    TPZManVector<int64_t,4> nodes(4);
    nodes[0]=0; nodes[1]=1; nodes[2]=2; nodes[3]=3;
    int64_t idx=0;
    gmesh->CreateGeoElement(EQuadrilateral, nodes, matId, idx);

    // lados (para BCs)
    idx=1; { TPZVec<long> L(2); L[0]=0; L[1]=1; new TPZGeoElRefPattern<pzgeom::TPZGeoLinear>(idx, L, -1, *gmesh); }
    idx=2; { TPZVec<long> L(2); L[0]=2; L[1]=3; new TPZGeoElRefPattern<pzgeom::TPZGeoLinear>(idx, L, -2, *gmesh); }

    gmesh->BuildConnectivity();

    for (int r=0;r<nRef;r++){
        const int nel = gmesh->NElements();
        TPZManVector<TPZGeoEl*> sub;
        for (int i=0;i<nel;i++) if (auto *eg = gmesh->ElementVec()[i]) eg->Divide(sub);
    }
    return gmesh;
}

// ------------------------------------------------------------
// Malha elástica H1 com memória (TPZMatElastic2DMem<TPZElasticMem>)
// ------------------------------------------------------------
TPZCompMesh * CompMeshElastic(TPZGeoMesh * gmesh, int porder, int matId)
{
    auto *cmesh = new TPZCompMesh(gmesh);
    cmesh->SetDimModel(2);
    cmesh->SetDefaultOrder(porder);
    cmesh->SetAllCreateFunctionsContinuousWithMem();

    auto *mat = new TPZMatElastic2DMem<TPZElasticMem>(matId);
    mat->SetId(1);      // (mantido)
    mat->SetId(matId);  // (mantido)
    REAL E=1., nu=0.;
    TPZElasticResponse ER; ER.SetEngineeringData(E, nu);
    mat->SetUpdateMem(true);
    mat->SetElasticResponse(ER);
    mat->SetUpdateMem(false);

    cmesh->InsertMaterialObject(mat);

    TPZFMatrix<STATE> val1 ( 2,2,0. );
    TPZVec<STATE> val2 ( 2,0. );

    // engasta inferior
    val2[0]=1.; val2[1]=1.;
    auto bcclamp = mat->CreateBC(mat,-1,3,val1,val2);
    // carga distribuída no topo (exemplo)
    val2[0]=0.; val2[1]=1.;
    auto bcload  = mat->CreateBC(mat,-2,1,val1,val2);

    cmesh->InsertMaterialObject(bcclamp);
    cmesh->InsertMaterialObject(bcload);
    cmesh->AutoBuild();
    return cmesh;
}

// ------------------------------------------------------------
// Malha para o problema de autovalor do kernel KL
// ------------------------------------------------------------
TPZCompMesh * BuildCompMeshKL(TPZGeoMesh *gmesh, int porder, int matId)
{
    auto *cmesh = new TPZCompMesh(gmesh);
    cmesh->SetDimModel(2);
    cmesh->SetDefaultOrder(porder);
    cmesh->SetAllCreateFunctionsContinuous();

    const REAL Lx=1., Ly=1.;
    auto *mat = new TPZMatKLKernel(matId, 2, Lx, Ly);
    mat->SetId(matId);

    // modo analítico só para referência visual (opcional)
    constexpr double A = 1.15021, k = 1.30654;
    mat->SetExact([A,k](const TPZVec<REAL>& x, STATE& u, TPZFMatrix<STATE>& du){
        const double xx=x[0], yy=x[1];
        u = (STATE)(A*std::cos(k*xx)*std::cos(k*yy));
        du.Resize(2,1);
        du(0,0) = (STATE)(-A*k*std::sin(k*xx)*std::cos(k*yy));
        du(1,0) = (STATE)(-A*k*std::cos(k*xx)*std::sin(k*yy));
    });

    cmesh->InsertMaterialObject(mat);
    cmesh->AutoBuild();
    cmesh->AdjustBoundaryElements();
    cmesh->CleanUpUnconnectedNodes();
    return cmesh;
}

// ------------------------------------------------------------
// ComputeField: grava H(x) nos IPs da malha target (autovetores em LINHAS)
// ------------------------------------------------------------
void ComputeField(TPZCompMesh* source,TPZCompMesh* target,
                  TPZFMatrix<REAL> colfield)
{
    // material com memória
    auto *targetmatwithmem =
    dynamic_cast<TPZMatWithMem<TPZElasticMem>*>(target->MaterialVec()[1]);
    if(!targetmatwithmem){ std::cout<<"material com memoria nao inicializado\n"; DebugStop(); }
    targetmatwithmem->SetUpdateMem(true);
    TPZAdmChunkVector<TPZElasticMem> &mem = *targetmatwithmem->GetMemory();

    source->LoadSolution(colfield);
    source->CleanUpUnconnectedNodes();
    source->InitializeBlock();
    source->ExpandSolution();
    source->LoadReferences();

    // --- parâmetros da distribuição (ESCOLHA UMA):
    // (a) Normal: H = mu + sigma * g
    const REAL mu    = 1.0;
    const REAL sigma = 0.2;
    const bool use_lognormal = false;

    // (b) Lognormal: H = exp(lambda + xi * g)  (apenas se use_lognormal=true)
    const REAL lambda_ln = std::log(mu) - 0.5*std::log(1.0 + (sigma*sigma)/(mu*mu)); // ex. coerente com (mu,sigma)
    const REAL xi_ln     = std::sqrt(std::log(1.0 + (sigma*sigma)/(mu*mu)));

    const int nels = target->NElements();
    for (int iel=0; iel<nels; iel++)
    {
        auto *targetcel = target->ElementVec()[iel];
        if(!targetcel) continue;
        auto *targetintel = dynamic_cast<TPZInterpolationSpace*>(targetcel);
        if(!targetintel) continue;

        auto *intelmat = dynamic_cast<TPZMatWithMem<TPZElasticMem>*>(targetintel->Material());
        if (intelmat != targetmatwithmem) continue;

        TPZMaterialDataT<STATE> datatarget,datasource;

        TPZIntPoints& intrule = targetintel->GetIntegrationRule();
        const int nip = intrule.NPoints();

        for (int ip=0; ip<nip; ip++)
        {
            REAL w;
            TPZManVector<REAL,3> intpt(3,0.);
            intrule.Point(ip, intpt, w);

            datatarget.intLocPtIndex = ip;
            datatarget.fNeedsSol = false;
            targetintel->InitMaterialData(datatarget);
            targetintel->ComputeRequiredData(datatarget, intpt);

            // acha, na malha source, o elemento e o ponto equivalente
            TPZManVector<REAL,3> qsisource(2,0.);
            int64_t elidsrc;
            TPZGeoEl *gelsource =source->Reference()->FindElement(datatarget.x, qsisource, elidsrc, source->Dimension());
            if (!gelsource) DebugStop();

            auto *celsource = gelsource->Reference();
            if (!celsource) DebugStop();
            auto *intelsource = dynamic_cast<TPZInterpolationSpace*>(celsource);
            if (!intelsource) DebugStop();

            datasource.intLocPtIndex = ip;
            datasource.fNeedsSol = true;
            intelsource->InitMaterialData(datasource);
            intelsource->ComputeRequiredData(datasource, qsisource);

            REAL H;
            if (!use_lognormal) {
                H = mu + sigma * datasource.sol[0][0];                 // NORMAL
            } else {
                H = std::exp(lambda_ln + xi_ln * datasource.sol[0][0]); // LOGNORMAL

            }
            if (H <= 0.) H = std::numeric_limits<REAL>::min();
            //std::cout<<"H = "<<H  <<std::endl;
            const int indextarget = datatarget.intGlobPtIndex;
            TPZElasticResponse ER;
            ER.SetEngineeringData(/*E=*/H, /*nu=*/0.0);
            mem[indextarget].m_ER = ER;
        }
    }
    targetmatwithmem->SetUpdateMem(false);
}

// ------------------------------------------------------------
// Resolve sistema elástico e retorna uy no ponto
// ------------------------------------------------------------
static REAL SolveOnce_Uy(TPZCompMesh* cmeshtarget, REAL xprobe, REAL yprobe)
{
    TPZLinearAnalysis an(cmeshtarget);
    TPZSkylineStructMatrix<STATE> str(cmeshtarget);
    an.SetStructuralMatrix(str);
    TPZStepSolver<REAL> direct; direct.SetDirect(ELDLt);
    an.SetSolver(direct);
    an.Assemble();
    an.Solve();
    return UyAtPoint(cmeshtarget, xprobe, yprobe);
}

// ------------------------------------------------------------
// MAIN Monte Carlo
// ------------------------------------------------------------
int main()
{
    const int   kMatId   = 1;
    const int   kPOrder  = 3;
    const int   kRef     = 3;
    const int   M        = 4;    // truncagem KL
    const int   Nsamp    = 30000;    // nº amostras MC
    const REAL  xprobe   = -0.5;  // ponto onde medimos uy
    const REAL  yprobe   =  0.5;

    // malhas
    TPZGeoMesh*  gmesh       = CreateGeoMeshMathematicaLike(kMatId, kRef);
    TPZCompMesh* cmeshKL     = BuildCompMeshKL(gmesh, kPOrder, kMatId);
    TPZCompMesh* cmeshtarget = CompMeshElastic( gmesh, kPOrder, kMatId);

    TPZFMatrix<REAL>  hhat;
    int runklsim=0;
    if(runklsim)
    {
        // autovalores/autovetores do kernel (Galerkin + Krylov)
        TPZEigenAnalysis an(cmeshKL,false);
        pzdoublestrmatriz<REAL> sm(cmeshKL);
        sm.SetCAssembly(pzdoublestrmatriz<REAL>::ECAssembly::Galerkin);
        an.SetStructuralMatrix(sm);

        TPZKrylovEigenSolver<STATE> esolver;
        esolver.SetAsGeneralised(true);
        esolver.SetEigenSorting(TPZEigenSort::AbsDescending);
        esolver.SetNEigenpairs(M);
        esolver.SetKrylovDim(cmeshKL->NEquations());
        esolver.SetTolerance(1e-10);
        an.SetSolver(esolver);

        an.Assemble();
        an.Solve();

        TPZFMatrix<CSTATE> evecs = an.Eigenvectors();//sao armazenados em colunas
        TPZVec<CSTATE>     evals = an.Eigenvalues();

        REAL cov=0.2;
        REAL mean=1;

        TPZFMatrix<REAL> THETA =ComputeTheta ( M,  mean,  cov, Nsamp);
        TPZFMatrix<REAL>  PHI=BuildPhiSqrtLambda(cmeshKL,evecs,evals,M);
        PHI.Multiply ( THETA, hhat );
        TPZBFileStream out;
        out.OpenWrite("hhat.bin");
        hhat.Write(out,0);
    }else{
        TPZBFileStream in;
        in.OpenRead("hhat.bin");
        hhat.Read(in,0);
    }

    std::ofstream out("mc_results.csv");
    out.setf(std::ios::fixed); out << std::setprecision(10);
    out << "sample,uy\n";
    double sum=0.0, sum2=0.0;
    int samp2=20;
    for (int s=0; s<samp2; ++s)
    {
        TPZFMatrix<REAL> colfield(hhat.Rows(),1);
        for (int k=0; k<hhat.Rows(); k++) colfield(k,0) = hhat(k,s);
        ComputeField(cmeshKL, cmeshtarget,colfield);
        const REAL uy = SolveOnce_Uy(cmeshtarget, xprobe, yprobe);
        out << s << "," << uy << "\n";
        sum  += uy;
        sum2 += (double)uy*(double)uy;
    }
    out.close();

    const double meanx = sum/samp2;
    const double var  = std::max(0.0, sum2/samp2 - meanx*meanx);
    const double stdv = std::sqrt(var);
    std::cout << "Monte Carlo: samp2="<<samp2
    << "  mean(uy)="<<meanx
    << "  std(uy)="<<stdv << std::endl;

    // VTK da última amostra
    PostElastic(cmeshtarget, "elastic_last.vtk", /*matid=*/1);

    delete cmeshtarget;
    delete cmeshKL;
    delete gmesh;
    return 0;
}
