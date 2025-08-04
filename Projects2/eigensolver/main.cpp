#include "TPZEigenAnalysis.h"
#include "TPZKrylovEigenSolver.h"
#include "TPZLapackEigenSolver.h"
#include "pzgeoelbc.h"
#include <pzgmesh.h> //for TPZGeoMesh
#include <pzcmesh.h> //for TPZCompMesh
#include "TPZLinearAnalysis.h"
//#include <catch2/catch.hpp>
#include "Poisson/TPZMatPoisson.h"
#include "pzmanvector.h"
#include <iostream>

// Função para criar uma malha 1D simples
TPZGeoMesh *Create1DGeoMesh(int nel, REAL L) {
        TPZGeoMesh *gmesh = new TPZGeoMesh();
        gmesh->NodeVec().Resize(nel+1);
        for(int i=0; i<=nel; i++) {
                gmesh->NodeVec()[i].SetNodeId(i);
                gmesh->NodeVec()[i].SetCoord(0, L*i/nel);
        }
        for(long i=0; i<nel; i++) {
                TPZVec<int64_t> indices(2);
                indices[0] = i;
                indices[1] = i+1;
                gmesh->CreateGeoElement(EOned, indices, 1, i);
        }
        gmesh->BuildConnectivity();
        return gmesh;
}

// Função para criar a malha computacional
TPZCompMesh *Create1DCompMesh(TPZGeoMesh *gmesh, int pOrder) {
        TPZCompMesh *cmesh = new TPZCompMesh(gmesh);
        int dim = 1;
        cmesh->SetDefaultOrder(pOrder);
        cmesh->SetDimModel(dim);

        // Material do tipo Laplace (autovalor)
        long matid = 1;
        TPZMatPoisson<REAL> *mat = new TPZMatPoisson<REAL>(matid, dim);
        //mat->SetSymmetric();
        cmesh->InsertMaterialObject(mat);

        // Condições de contorno homogêneas
        TPZFMatrix<STATE> val1(1,1,0.); // matriz vazia
        TPZManVector<STATE,1> val2(1,0.);

        int bc_left = -1, bc_right = -2;
        auto bcleft = mat->CreateBC(mat, bc_left, 0, val1, val2);
        auto bcright = mat->CreateBC(mat, bc_right, 0, val1, val2);

        //TPZBndCond *bcleft = mat->CreateBC ( mat,bc_left,0,val1,val2 ); //clamped line
        //TPZBndCond * bcright = mat->CreateBC(mat, bc_right, 0, val1, val2);

        cmesh->InsertMaterialObject(bcleft);
        cmesh->InsertMaterialObject(bcright);

        // Associa os BCs geométricos
        gmesh->Element(0)->SetMaterialId(bc_left);
        gmesh->Element(gmesh->NElements()-1)->SetMaterialId(bc_right);

        cmesh->AutoBuild();
        return cmesh;
}

int main() {
        //Parâmetros do problema
        int nel = 10;         // número de elementos
        REAL L = 1.0;         // comprimento
        int pOrder = 2;       // ordem do polinômio

        // Monta as malhas
        TPZGeoMesh *gmesh = Create1DGeoMesh(nel, L);
        TPZCompMesh *cmesh = Create1DCompMesh(gmesh, pOrder);

        // Análise de autovalores
        TPZEigenAnalysis analysis(cmesh);

        // Solver de autovalores (LAPACK)
        //TPZLapackEigenSolver<CSTATE> solver;
        TPZKrylovEigenSolver<CSTATE> solver;
        solver.SetNEigenpairs(3); // calcular 3 autovalores

        //analysis.SetSolver(solver);

        // Montar as matrizes
        analysis.Assemble();

        // Resolver
        analysis.Solve();

        // Pega os autovalores/vetores
        //TPZVec<CSTATE> eigval = analysis.GetEigenvalues();
        //TPZFMatrix<CSTATE> eigvec = analysis.GetEigenvectors();

        std::cout << "Autovalores:" << std::endl;
        //for(int i=0; i<eigval.size(); i++) {
                //std::cout << eigval[i] << std::endl;
        //}

        std::cout <<"HELLO WORLD"<<std::endl;

        return 0;
}
