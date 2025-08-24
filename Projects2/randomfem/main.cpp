// main.cpp — Passo 1: KL + geração de campo lognormal E(x)
// Compila com NeoPZ (usa TPZEigenAnalysis/TPZKrylovEigenSolver/TPZMatKLKernel)

#include "pzgmesh.h"
#include "pzcmesh.h"
#include "pzgeoquad.h"
#include "tpzgeoelrefpattern.h"
#include "pzgeoelbc.h"
#include "TPZGeoLinear.h"
#include "TPZVTKGeoMesh.h"
#include "pzintel.h"
#include "pzstack.h"

#include "TPZEigenAnalysis.h"
#include "TPZKrylovEigenSolver.h"
#include "pzdoublestrmatriz.h"
#include "pzskylstrmatrix.h"

#include "TPZMatKLKernel.h"
#include "Elasticity/TPZMatElastic2DMem.h"
#include "Elasticity/TPZElasticMem.h"
#include <random>
#include <cmath>
#include <algorithm>
#include <vector>
#include <string>

// ---------------- geometria simples (quadrado [-0.5,0.5]^2), BC geométricas só para exemplo VTK ----------------
static TPZAutoPointer<TPZGeoMesh>
CreateGeoMeshMathematicaLike(int matId, int nRef)
{
        TPZAutoPointer<TPZGeoMesh> gmesh = new TPZGeoMesh;
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

        // BC geométrica (apenas para ter fronteiras nomeadas no VTK)
        idx=1;
        TPZVec<long> TopoLine(2); TopoLine[0]=0; TopoLine[1]=1;
        new TPZGeoElRefPattern<pzgeom::TPZGeoLinear>(idx, TopoLine, -1, *gmesh);

        idx=2;
        TopoLine[0]=2; TopoLine[1]=3;
        new TPZGeoElRefPattern<pzgeom::TPZGeoLinear>(idx, TopoLine, -2, *gmesh);

        gmesh->BuildConnectivity();

        for (int r=0;r<nRef;r++){
                const int nel = gmesh->NElements();
                TPZManVector<TPZGeoEl*> sub;
                for (int i=0;i<nel;i++) if (auto *eg = gmesh->ElementVec()[i]) eg->Divide(sub);
        }
        return gmesh;
}

// ---------------- CompMesh H1 para o problema de autovalor do kernel (KL) ----------------
static TPZCompMesh*
CompMeshElastic(TPZAutoPointer<TPZGeoMesh> gmesh, int porder, int matId)
{
        auto *cmesh = new TPZCompMesh(gmesh);
        unsigned int dim  = 2;
        const std::string name ( "ElastoPlastic COMP MESH Footing Problem " );
        cmesh->SetName ( name );
        cmesh->SetDimModel(dim);
        cmesh->SetDefaultOrder(porder);
        cmesh->SetAllCreateFunctionsContinuousWithMem();

        auto *mat = new TPZMatElastic2DMem<TPZElasticMem>(matId);
        mat->SetId ( 1 );
        mat->SetId(matId);
        REAL E=1.;
        REAL nu=0.;
        mat->SetElasticityFallback(E,  nu);





        // Elastic predictor
        TPZElasticResponse ER;
        //TPZPlasticStepPV<TPZYCMohrCoulombPV, TPZElasticResponse> LEMC;
        //TPZPlasticStepVoigt<TPZMohrCoulombVoigt, TPZElasticResponse> LEMC;
        ER.SetEngineeringData( E, nu );
        mat->SetElasticResponse ( ER );


        cmesh->InsertMaterialObject(mat);

        TPZFMatrix<STATE> val1 ( 2,2,0. );
        TPZVec<STATE> val2 ( 2,0. );

        val2[0]=1.;
        val2[1]=1.;
        auto bcclamp = mat->CreateBC ( mat,-1,3,val1,val2 ); //clamped line restrictions

        val2[0]=0.;
        val2[1]=1.;
        auto bcload = mat->CreateBC ( mat,-2,1,val1,val2 ); //distributed load superior line


        cmesh->InsertMaterialObject(bcclamp);
        cmesh->InsertMaterialObject(bcload);
        cmesh->AutoBuild();
        return cmesh;
}


// ---------------- CompMesh H1 para o problema de autovalor do kernel (KL) ----------------
static TPZCompMesh*
BuildCompMeshKL(TPZAutoPointer<TPZGeoMesh> gmesh, int porder, int matId)
{
        auto *cmesh = new TPZCompMesh(gmesh);
        cmesh->SetDimModel(2);
        cmesh->SetDefaultOrder(porder);
        cmesh->SetAllCreateFunctionsContinuous();

        const REAL Lx=1., Ly=0.1;
        auto KernelFn = [Lx,Ly](const TPZVec<REAL>& x, const TPZVec<REAL>& y)->STATE{
                return std::exp(-std::fabs(x[0]-y[0])/Lx - std::fabs(x[1]-y[1])/Ly);
        };
        auto *mat = new TPZMatKLKernel(matId, 2, KernelFn);
        mat->SetId(matId);

        // (opcional) modo analítico para visualização do modo 0
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


TPZFMatrix<REAL> CreateLogNormalRandomField ( TPZFMatrix<REAL> PHI, REAL mean, REAL cov,int samples)
{

        int M = PHI.Cols();

        std::normal_distribution<double> distribution ( 0., 1. );

        std::random_device rd{};

        std::mt19937 generator{ rd() };


        TPZFMatrix<REAL>  THETA ( M, samples, 0. );

        for ( int isample = 0; isample < samples; isample++ )
        {
                for ( int irdvar = 0; irdvar < M; irdvar++ )
                {

                        REAL xic = distribution ( generator );
                        THETA ( irdvar,isample ) = xic;
                }
        }


        TPZFMatrix<REAL>  hhat;

        PHI.Multiply ( THETA, hhat );

        REAL sdev = cov * mean;
        REAL xi = sqrt ( log ( 1 + pow ( ( sdev / mean ),2 ) ) );
        REAL lambda = log ( mean ) - xi * xi / 2.;
        for ( int i = 0; i < hhat.Rows(); i++ )
        {
                for ( int j = 0; j < hhat.Cols(); j++ )
                {
                        REAL temp =  hhat ( i,j );
                        hhat ( i,j ) = exp ( lambda + xi * temp );
                }
        }


        return hhat;
}

void TransferSolutionFrom(TPZCompMesh* cmeshElastic, TPZCompMesh* cmeshKL, int matidElastic = 1)
{
        // 1) material elástico com memória no TARGET
        auto* matEl = dynamic_cast<TPZMatElastic2DMem<TPZElasticMem>*>(cmeshElastic->FindMaterial(matidElastic));
        if (!matEl) DebugStop();
        auto* wm = static_cast<TPZMatWithMem<TPZElasticMem>*>(matEl);

        // garante alocação da memória e habilita escrita
        wm->SetUpdateMem(true);
        // nome pode variar por branch: InitializeElementMemory / FillMemories / InitializeMemory
        //wm->InitializeElementMemory();

        // 2) preparar a KL para busca de elemento contendo x
        TPZGeoMesh* gmeshKL = cmeshKL->Reference();
        const int dim = gmeshKL->Dimension();
        gmeshKL->ResetReference();
        cmeshKL->LoadReferences();

        // 3) varrer elementos do TARGET
        const int64_t nels = cmeshElastic->NElements();
        for (int64_t iel = 0; iel < nels; ++iel)
        {
                auto* celTar = cmeshElastic->ElementVec()[iel];
                auto* intelTar = dynamic_cast<TPZInterpolationSpace*>(celTar);
                if (!intelTar) continue;

                TPZMaterial* base = intelTar->Material();
                if (!base || base->Id() != matidElastic) continue;

                TPZMaterialDataT<STATE> dataTar;
                intelTar->InitMaterialData(dataTar);

                const TPZIntPoints& rule = intelTar->GetIntegrationRule();
                const int nip = rule.NPoints();

                for (int ip = 0; ip < nip; ++ip)
                {
                        // (a) dados no ALVO (para obter x e o índice global do IP na memória do elástico)
                        REAL w; TPZManVector<REAL,3> qsiTar(3,0.);
                        rule.Point(ip, qsiTar, w);

                        dataTar.intLocPtIndex = ip;
                        dataTar.fNeedsSol = false;               // só precisamos de x e índices aqui
                        intelTar->ComputeRequiredData(dataTar, qsiTar);

                        const int64_t idMem = dataTar.intGlobPtIndex;  // índice global do IP do ALVO

                        // (b) localizar ponto na KL e avaliar E(x)
                        TPZManVector<REAL,3> qsiKL(3,0.);
                        long elidKL;
                        TPZGeoEl* gelKL = gmeshKL->FindElement(dataTar.x, qsiKL, elidKL, dim);
                        if (!gelKL) DebugStop();

                        auto* celKL  = gelKL->Reference();
                        auto* intelKL = dynamic_cast<TPZInterpolationSpace*>(celKL);
                        if (!intelKL) DebugStop();

                        TPZMaterialDataT<STATE> dataKL;
                        intelKL->InitMaterialData(dataKL);
                        dataKL.fNeedsSol = true;                 // aqui sim precisamos da solução
                        intelKL->ComputeRequiredData(dataKL, qsiKL);

                        const REAL E = (dataKL.sol.size() && dataKL.sol[0].size()) ? dataKL.sol[0][0] : 0.0;

                        // (c) escrever na memória do ALVO
                        TPZElasticMem& m = wm->MemItem(idMem);
                        const REAL nu_keep = m.m_ER.Poisson();   // mantém ν que já estava (ou defina o seu)
                        m.m_ER.SetEngineeringData(E, nu_keep);
                }
        }

        wm->SetUpdateMem(false);
}
// ElasticityTools.h
#include "pzpostprocanalysis.h"
void PostElastic(TPZCompMesh* cmesh,const std::string& vtkfile,int matid = 1);
void CreatePostProcessingMesh(TPZCompMesh* cmesh,TPZPostProcAnalysis* pproc,int matid = 1);
#include "pzstepsolver.h"
void SolveElastic(TPZCompMesh * source)
{
        int matid=1,porder=2,ref=3;
        TPZAutoPointer<TPZGeoMesh> gmesh  = CreateGeoMeshMathematicaLike(matid, ref);
        TPZCompMesh* cmeshelastic = CompMeshElastic(gmesh, porder,  matid);

        TransferSolutionFrom(cmeshelastic,  source);

        auto anal= TPZLinearAnalysis(cmeshelastic);
        TPZSkylineStructMatrix<STATE> strmat(cmeshelastic);
        anal.SetStructuralMatrix(strmat);
        auto direct = new TPZStepSolver<REAL>;
        direct->SetDirect ( ELDLt );
        anal.SetSolver ( *direct );
        anal.Assemble();

        anal.Rhs().Print("RHS");


        // auto base = anal.Solver();
        // if (auto step = dynamic_cast<TPZStepSolver<STATE>*>(base)) {
        //         TPZAutoPointer<TPZMatrix<STATE>> K = step->Matrix(); // nem todas as versões têm isso
        //         if (K) K->Print("K = ", std::cout);
        // }
        //anal.Matrix();
        anal.Solve();

        anal.Solution().Print("sol");

        PostElastic(cmeshelastic, "Elasticity.vtk", /*matid=*/1);

        //auto K = anal.fSolver;
        //auto matrix=anal.Solver()->Matrix();

        // (opcional) VTK do modo 0, parte real
        {
               //TPZStack<std::string> scalars, vectors;
                //vectors.Push("Displacement");
                //scalars.Push("Young");
                //scalars.Push("Poisson");
                //anal.DefineGraphMesh(2, scalars, vectors, "Elasticity.vtk");
              // anal.PostProcess(0);
        }

}

int main()
{



        // ---- parâmetros gerais ----
        static const int      kMatId   = 1;
        static const int      kPOrder  = 2;
        static const int      kRef     = 3;       // níveis de refino geométrico


        // ---- malha + CMesh do KL ----
        TPZAutoPointer<TPZGeoMesh> gmesh  = CreateGeoMeshMathematicaLike(kMatId, kRef);
        TPZCompMesh* cmeshKL              = BuildCompMeshKL(gmesh, kPOrder, kMatId);

        // ---- análise de autovalor (Galerkin + Krylov) ----
        TPZEigenAnalysis an(cmeshKL);

        pzdoublestrmatriz<STATE> sm(cmeshKL);
        sm.SetCAssembly(pzdoublestrmatriz<STATE>::ECAssembly::Galerkin);
        an.SetStructuralMatrix(sm);

        const int ndof = cmeshKL->NEquations();

        TPZKrylovEigenSolver<STATE> solver;
        solver.SetAsGeneralised(true);
        solver.SetEigenSorting(TPZEigenSort::AbsDescending);
        solver.SetNEigenpairs(ndof);
        solver.SetKrylovDim(ndof);
        solver.SetTolerance(1e-10);
        an.SetSolver(solver);

        an.Assemble();
        an.Solve();

        // (opcional) VTK do modo 0, parte real
        {
                TPZStack<std::string> scalars, vectors;
                scalars.Push("Solution");
                scalars.Push("ExactSolution");
                vectors.Push("Gradient");
                vectors.Push("ExactGradient");
                an.DefineGraphMesh(2, scalars, vectors, "kl_mode0.vtk");
                an.PostProcessMode(0, 0, scalars, vectors, TPZEigenAnalysis::EEigPart::Real, true);
        }

         const int M = 20; // por exemplo

        an.CheckL2Norms(M,std::cout);

        an.BuildPhiSqrtLambdaNodal(M);

        TPZFMatrix<STATE>  sqrtvalvec=an.Solution();

        REAL mean=1.;

        REAL cov=0.3;

        int samples=3;

        TPZFMatrix<STATE> field = CreateLogNormalRandomField(sqrtvalvec,mean,cov,samples);

        for(int i=0;i<2;i++)
        {
                TPZFMatrix<REAL> colmat(field.Rows(),1);
                for(int j=0;j<field.Rows();j++)colmat(j,0)=field(j,i);
                cmeshKL->LoadSolution(colmat);
                SolveElastic(cmeshKL);

                TPZStack<std::string> scalars, vectors;
                scalars.Push("Solution");
                an.DefineGraphMesh(2, scalars, vectors, "campos.vtk");
                an.PostProcess(0);
                //an.PostProcessMode(0, 0, scalars, vectors, TPZEigenAnalysis::EEigPart::Real, true);
        }

        delete cmeshKL;
        return 0;
}




// nomes que o seu material entende em VariableIndex(...)
void PostProcessVariables(TPZStack<std::string>& scal, TPZStack<std::string>& vec)
{
                // escalares
        scal.Push("Young");
        scal.Push("Poisson");
        scal.Push("POrder");
                // vetores (3 componentes cada)
        vec.Push("Displacement"); // (ux,uy,0)
        vec.Push("Strain");       // (exx, eyy, gxy)
        vec.Push("Stress");       // (sxx, syy, sxy)
}

        // cria/atualiza a malha de pós-processo e transfere a solução
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
                // copia a última solução do cmesh para a malha de pós-processo
        pproc->TransferSolution();
}

        // cria o VTK (dim=2) com as variáveis acima
void PostElastic(TPZCompMesh* cmesh,const std::string& vtkfile,int matid)
{
        TPZPostProcAnalysis pproc;
        CreatePostProcessingMesh(cmesh, &pproc, matid);

        TPZStack<std::string> scal, vec;
        PostProcessVariables(scal, vec);

        pproc.DefineGraphMesh(/*dim=*/2, scal, vec, vtkfile);
        pproc.PostProcess(0);
}

