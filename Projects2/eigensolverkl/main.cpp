// main.cpp — KL em [-0.5,0.5]^2 (3x3 quads) com Nyström explícito + TPZKrylovEigenSolver

#include "pzgmesh.h"
#include "pzcmesh.h"
#include "pzgeoquad.h"
#include "pzintel.h"

#include "TPZEigenAnalysis.h"
#include "TPZKrylovEigenSolver.h"
#include "TPZLapackEigenSolver.h"
#include "pzskylstrmatrix.h"

#include "TPZMatKLCov2D.h"
#include "TPZVTKGeoMesh.h"

#include <iostream>
#include <algorithm>
#include <complex>
#include <fstream>

// ------------------ malha [-0.5,0.5]^2 com nx×ny QUADs ------------------
#include "pzgmesh.h"
#include "pzgeoquad.h"
#include "TPZVTKGeoMesh.h"
#include "pzdoublestrmatriz.h"
#include "TPZMatKLKernel.h"
// malha [-0.5,0.5]x[-0.5,0.5], nós e conectividades iguais às do Mathematica
// matId: id do material geométrico a ser atribuído aos elementos
TPZGeoMesh* CreateGeoMeshMathematicaLike(int matId)
{
        auto gmesh = new TPZGeoMesh;
        gmesh->SetDimension(2);

        // ---- nós (mesma ordem que o Mathematica) ----
        // 1..16  ->  0..15 aqui
        const double m = 1.0/6.0;
        const double X[16][2] = {
                {-0.5, -0.5},    {-0.5, -m},      {-0.5,  m},      {-0.5,  0.5},
                {-m , -0.5},     {-m , -m},       {-m ,  m},       {-m ,   0.5},
                { m, -0.5},      { m, -m},        { m,  m},        { m,   0.5},
                {0.5, -0.5},     {0.5, -m},       {0.5,  m},       {0.5,  0.5}
        };

        gmesh->NodeVec().Resize(16);
        for (int i = 0; i < 16; i++) {
                TPZManVector<REAL,3> coord(3, 0.);
                coord[0] = X[i][0];
                coord[1] = X[i][1];
                gmesh->NodeVec()[i].Initialize(coord, *gmesh);
        }

        // ---- conectividade (1-based do Mathematica -> 0-based aqui) ----
        const int conn[9][4] = {
                { 1, 5, 6, 2}, { 2, 6, 7, 3}, { 3, 7, 8, 4},
                { 5, 9,10, 6}, { 6,10,11, 7}, { 7,11,12, 8},
                { 9,13,14,10}, {10,14,15,11}, {11,15,16,12}
        };

        for (int e = 0; e < 9; e++) {
                TPZManVector<int64_t,4> nodes(4);
                nodes[0] = conn[e][0] - 1;
                nodes[1] = conn[e][1] - 1;
                nodes[2] = conn[e][2] - 1;
                nodes[3] = conn[e][3] - 1;

                int64_t index;
                gmesh->CreateGeoElement(EQuadrilateral, nodes, matId, index);
        }

        gmesh->BuildConnectivity();

        // (opcional) exporta para VTK para você conferir
        std::ofstream vtk("gmesh_from_mathematica.vtk");
        TPZVTKGeoMesh::PrintGMeshVTK(gmesh, vtk, false);

        return gmesh;
}


// ------------------ CompMesh H1 ------------------
static TPZCompMesh* MakeCompMesh(TPZAutoPointer<TPZGeoMesh> gmesh, int porder, int matId)
{
        auto *cmesh = new TPZCompMesh(gmesh);
        cmesh->SetDimModel(gmesh->Dimension());
        cmesh->SetDefaultOrder(porder);
        cmesh->SetAllCreateFunctionsContinuous();

        auto *mat = new TPZMatKLCov2D(matId, gmesh->Dimension());
        // parâmetros do kernel (Lx, Ly, sigma^2)
        mat->SetKernelParams(1.0, 1.0, 1.0);

        cmesh->InsertMaterialObject(mat);
        cmesh->AutoBuild();
        cmesh->AdjustBoundaryElements();
        cmesh->CleanUpUnconnectedNodes();
        return cmesh;
}
// ------------------ CompMesh H1 ------------------
static TPZCompMesh* CompMesh(TPZAutoPointer<TPZGeoMesh> gmesh, int porder, int matId)
{
        auto *cmesh = new TPZCompMesh(gmesh);
        cmesh->SetDimModel(gmesh->Dimension());
        cmesh->SetDefaultOrder(porder);
        cmesh->SetAllCreateFunctionsContinuous();
        REAL Lx=1.;
        REAL Ly=1.;
        auto KernelFn=[Lx,Ly](const TPZVec<REAL> &x,const TPZVec<REAL>&y)->STATE{
                const REAL dx=x[0]-y[0];
                const REAL dy=x[1]-y[1];
                REAL val = exp ( -fabs ( dx ) / ( Lx )- fabs (dy ) / ( Ly ) );
                return val;
        };


        //auto * mat = new TPZElasticity2D ( 1,21000000.,0.3,0.,0. ); //selfweigth
        auto *mat = new TPZMatKLKernel(matId,2, KernelFn);

        const bool isGEV =  (dynamic_cast<TPZMatGeneralisedEigenVal*>(mat)!=nullptr );

        std::cout <<" isGEV = "<<isGEV<<"\n";
        mat->SetId(matId);

        cmesh->InsertMaterialObject(mat);
        cmesh->AutoBuild();
        cmesh->AdjustBoundaryElements();
        cmesh->CleanUpUnconnectedNodes();
        return cmesh;
}

// ------------------ main ------------------
int main()
{
        const int   matId  = 1;
        const int   porder = 2;          // MeshOrder -> 1 (como no Mathematica)

        // 1) Malhas
        auto gmesh = CreateGeoMeshMathematicaLike(matId);
        TPZCompMesh *cmesh = CompMesh(gmesh, porder, matId);


        TPZEigenAnalysis an(cmesh);

        pzdoublestrmatriz<STATE> sm(cmesh);


        an.SetStructuralMatrix(sm);

        const int nact = cmesh->NEquations();


        //int solver= enum {ELAPCK=1,EKRILOV=2};
        //TPZKrylovEigenSolver<STATE> solver;
        TPZLapackEigenSolver<STATE> solver;
       // auto solver= new TPZLapackEigenSolver<STATE>;
        solver.SetAsGeneralised(true);
        solver.SetNEigenpairs(nact);
        //solver.SetKrylovDim(nact);
        //solver.SetTolerance(1e-12);
        solver.SetEigenSorting(TPZEigenSort::AbsDescending); // maiores autovalores primeiro

        an.SetSolver(solver);

        // for(auto &it :cmesh->MaterialVec()){
        //         TPZMaterial *m=it.second;
        //         const bool isBC =  (dynamic_cast<TPZBndCond*>(m)!=nullptr );
        //         const bool isGEV =  (dynamic_cast<TPZMatGeneralisedEigenVal*>(m)!=nullptr );
        //         std::cout << " id = "<<it.first << " type = "<<typeid(*m).name()<<" isBC "<<isBC <<" isGEV = "<<isGEV<<"\n";
        // }

        //std::cout << "Montando matrizes" << std::endl;
        an.Assemble();

        //solver.SetAsGeneralised(false);
        //auto &es = an.EigenSolver<STATE>();

        //solver.SetAsGeneralised(false);
        //auto A=es.MatrixA();
        //A->Print("MatrixA");
        //auto B=es.MatrixB();
        //B->Print("MatrixB");
        //
        // TPZFMatrix<REAL> invB,solM;
        // B->Inverse(invB,ELDLt);
        //
        // invB.Multiply(*A,solM);
        //
        // TPZAutoPointer<TPZMatrix<STATE>> Ap = new TPZFMatrix<STATE>(solM);
        //
        // solver.SetMatrixA(Ap);
        // std::cout << "aqui" << std::endl;
        // TPZFMatrix<CSTATE> vecs;
        // TPZVec<CSTATE> vals;
        // solver.SolveEigenProblem(vals,vecs);
        //
        //
        // for(int i=0;i<10;i++)std::cout << vals[i].real()<<std::endl;

         an.Solve();
         TPZVec<CSTATE> vals = an.GetEigenvalues();
         for(int i=0;i<10;i++)std::cout << vals[i].real()<<std::endl;

        // return 0;
        // if(true)
        // {
        //
        //
        // //cmesh->StructMatrix()->EquationFilter().Reset();
        //
        // // 2) B e C por Nyström explícito
        // TPZFMatrix<STATE> B, C;
        // const int qorder = 10; // regra "segura"
        //
        // TPZMatKLCov2D::BuildB_Nystrom(*cmesh, qorder, B);
        //
        // // kernel exponencial separável (Lx=Ly=1, sigma^2=1)
        // auto ker = [](const TPZVec<REAL>& x, const TPZVec<REAL>& y)->STATE {
        //         return TPZMatKLCov2D::ExpKernel(x,y,1.0,1.0,1.0);
        // };
        // TPZMatKLCov2D::BuildC_Nystrom(*cmesh, qorder, ker, C);
        //
        //
        // const int nact = cmesh->NEquations();
        //
        // std::cout << "B.Rows()" <<B.Rows()<<std::endl;
        // std::cout << "nact" <<nact<<std::endl;
        //
        // TPZKrylovEigenSolver<STATE> solver;
        // //solver.SetAsGeneralised(true);
        // solver.SetNEigenpairs(nact);
        // solver.SetKrylovDim(nact);
        // solver.SetTolerance(1e-12);
        // solver.SetEigenSorting(TPZEigenSort::AbsDescending); // maiores autovalores primeiro
        //
        //
        // TPZFMatrix<REAL> invB,solM;
        // B.Inverse(invB,ELDLt);
        //
        // invB.Multiply(C,solM);
        //
        // TPZAutoPointer<TPZMatrix<STATE>> Ap = new TPZFMatrix<STATE>(solM);
        //
        // solver.SetMatrixA(Ap);
        // std::cout << "aqui" << std::endl;
        // TPZFMatrix<CSTATE> vecs;
        // TPZVec<CSTATE> vals;
        // solver.SolveEigenProblem(vals,vecs);
        //
        //
        // for(int i=0;i<10;i++)std::cout << vals[i].real()<<std::endl;
        //
        //
        // }
        //


        return 0;
}
