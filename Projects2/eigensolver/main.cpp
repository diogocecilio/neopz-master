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
#include "TPZMatGeneralisedEigenVal.h"
#include "TPZVTKGeoMesh.h"
#include "tpzgeoelrefpattern.h"
#include "Elasticity/TPZElasticity2D.h"
#include "Elasticity/TPZElasticity2DGenEVP.h"
#include <pzskylstrmatrix.h>
#include "pzpostprocanalysis.h"
#include "Plasticity/TPZElasticResponse.h"
// Função para criar uma malha 1D simples
using namespace std;
TPZGeoMesh *CreateGeoMeshBending(int nref)
{
        // REAL co[4][2] = {{0.,0.},{5,0},{5,0.5},{0,0.5}};
        REAL co[12][2] = {{0., 0.}, {0., 0.5}, {1., 0.}, {1., 0.5}, {2., 0.}, {2., 0.5}, {3.,
                0.}, {3., 0.5}, {4., 0.}, {4., 0.5}, {5., 0.}, {5., 0.5}};
                long indices[5][4] = {{0,2,3,1},{2,4,5,3},{4,6,7,5},{6,8,9,7},{8,10,11,9}};
                TPZGeoEl *elvec[5];
                TPZGeoMesh *gmesh = new TPZGeoMesh();
                gmesh->SetDimension ( 2 );
                long nnode = 12;
                long nod;
                for ( nod=0; nod<nnode; nod++ )
                {
                        long nodind = gmesh->NodeVec().AllocateNewElement();
                        TPZVec<REAL> coord ( 2 );
                        coord[0] = co[nod][0];
                        coord[1] = co[nod][1];
                        gmesh->NodeVec() [nodind] = TPZGeoNode ( nod,coord,*gmesh );
                }

                long el;
                long nelem = 5;
                for ( el=0; el<nelem; el++ )
                {
                        TPZVec<long> nodind ( 4 );
                        for ( nod=0; nod<4; nod++ ) nodind[nod]=indices[el][nod];
                        long index;
                        elvec[el] = gmesh->CreateGeoElement ( EQuadrilateral,nodind,1,index);
                }

                TPZVec <long> TopoLine ( 2 );

                //long index=2;
                TopoLine[0] = 0;
                TopoLine[1] = 1;
                new TPZGeoElRefPattern< pzgeom::TPZGeoLinear> ( 5, TopoLine, - 1, *gmesh );//clamped in right side

                TPZVec <long> node ( 1 );
                node[0]=0;
                new TPZGeoElRefPattern< pzgeom::TPZGeoPoint> ( 6, node, - 2, *gmesh );//load node


                gmesh->BuildConnectivity();


                cout << "c" << endl;
                for ( int d = 0; d<nref; d++ )
                {
                        int nel = gmesh->NElements();
                        TPZManVector<TPZGeoEl *> subels;
                        for ( int iel = 0; iel<nel; iel++ )
                        {
                                TPZGeoEl *gel = gmesh->ElementVec() [iel];
                                gel->Divide ( subels );
                        }
                }

                //gmesh->Print(cout);
                std::ofstream files ( "teste-mesh.vtk" );
                TPZVTKGeoMesh::PrintGMeshVTK ( gmesh,files,false );
                cout << "d" << endl;
                return gmesh;

}
TPZCompMesh *CreateMeshBending ( TPZGeoMesh *gmesh,int poder )
{
        TPZCompMesh *cmesh = new TPZCompMesh ( gmesh );
        cmesh->SetDefaultOrder ( poder);

        //TPZElasticityMaterial(int id, REAL E, REAL nu, REAL fx, REAL fy, int planestress = 1);

        //auto * mat = new TPZElasticity2D ( 1,100.,0.,0.,0. ); //selfweigth
        auto * mat = new TPZElasticity2DGenEVP ( 1,100.,0.,0.,0., 100.,1.); //selfweigthTPZElasticity2DGenEVP

        cmesh->SetDimModel ( 1 );

        // TPZFMatrix<REAL> val1(2,2,0.),val2(2,1,0.);
        TPZFMatrix<STATE> val1 ( 2,2,0. );
        TPZVec<STATE> val2 ( 2,0. );
        //TPZMaterial *bcload,*bcclamp,*bcnode;

        val2[0]=0.;
        val2[1]=0.;
        auto bcclamp = mat->CreateBC ( mat,-1,0,val1,val2 ); //clamped line restrictions

        val2[0]=0.;
        val2[1]=0.;
        auto bcnode = mat->CreateBC ( mat,-2,0,val1,val2 ); //bottomrigth node restrictions


        cmesh->InsertMaterialObject ( mat );
        cmesh->InsertMaterialObject ( bcclamp );
        cmesh->InsertMaterialObject ( bcnode );


        cmesh->SetAllCreateFunctionsContinuous();

        cmesh->AutoBuild();
        cout << "e" << endl;
        cmesh->AdjustBoundaryElements();
        cmesh->CleanUpUnconnectedNodes();

        return cmesh;
}




int main() {

        // Monta as malhas
        int ref=0;
        int porder=1;
        auto gmesh = CreateGeoMeshBending(ref);
        auto cmesh = CreateMeshBending(gmesh,porder);

        // Análise de autovalores
        TPZEigenAnalysis analysis(cmesh);

        // criar struct matrix
        TPZSkylineStructMatrix<STATE> strmat(cmesh);
        strmat.SetNumThreads(0); // ou outro número


        analysis.SetStructuralMatrix(strmat);

        analysis.StructMatrix()->EquationFilter().Reset();
        // agora sim pode chamar
        const int nact = analysis.StructMatrix()->EquationFilter().NActiveEquations();

        cout << "numero de equações ativas = " << nact << endl;

        //quantos queremos (ex.: 8). Se não há CC, existem ~3 modos rígidos ~0.
        int nev = std::min(90, nact-1);   // não peça mais que n-1

        TPZKrylovEigenSolver<STATE> solver;
        solver.SetNEigenpairs(nev);
        // regra segura: KDim >= 2*nev + 10 e < nact
        int kdim = std::min(nact-1, std::max(30, 2*nev + 10));
        solver.SetKrylovDim(kdim);

        solver.SetTolerance(1e-10);
        solver.SetAsGeneralised(true);

        solver.SetEigenSorting(TPZEigenSort::AbsAscending);

        analysis.SetSolver(solver);
        analysis.Assemble();
        analysis.Solve();

        // 1) Autovalores/autovetores
        auto vals = analysis.GetEigenvalues();
        auto eig = analysis.GetEigenvectors();

        std::cout << "nev=" << vals.size()
        << "  eig(rows,cols)=(" << eig.Rows() << "," << eig.Cols() << ")\n";

        // 2) escolha do modo, sem estourar índice
        int mode = 0; // 0=1º, 1=2º, ...
        if (eig.Cols() == 0) { std::cerr << "Sem autovetores!\n"; return 0; }
        if (mode >= eig.Cols()) mode = eig.Cols()-1; // garante

        cout << "numero de autovetores "<<vals.size() << endl;
        for(int i=0;i<vals.size();i++)cout<< vals[i].real() <<endl;
        // 3) constrói solução com o TAMANHO CERTO (linhas = eig.Rows)
        TPZFMatrix<STATE> solRed(eig.Rows(), 1, 0.);
        for (int r = 0; r < eig.Rows(); r++) {
                solRed(r,0) = (STATE)eig(r,mode).real(); // ou imag()/abs()
        }

        // 4) espalha p/ o tamanho total se houver filtro de equações
        const int neqFull = analysis.Mesh()->NEquations();
        TPZFMatrix<STATE> solFull;
        if (eig.Rows() != neqFull) {
                // há redução: faz Scatter para o vetor completo
                analysis.StructMatrix()->EquationFilter().Scatter(solRed, solFull);
        } else {
                solFull = solRed;
        }

        // 5) carrega na análise e pós-processa
        analysis.Solution().Redim(solFull.Rows(), 1);
        analysis.Solution() = solFull;
        analysis.LoadSolution();

        TPZStack<std::string> scal, vec;
        scal.Push("EVP_U");           // ou "U"/"Ux"/"Uy" (conforme seu material)
        vec.Push("displacement");     // só se você implementou

        analysis.DefineGraphMesh(2, scal, vec, "ref3mode_0.vtk");
        analysis.PostProcess(0);

        return 0;
}
// int main() {
//
//         // Monta as malhas
//         TPZGeoMesh *gmesh = CreateGeoMeshBending();
//         TPZCompMesh *cmesh = CreateMeshBending(gmesh);
//
//         // Análise de autovalores
//         TPZEigenAnalysis analysis(cmesh);
//
//         // criar struct matrix
//         TPZSkylineStructMatrix<STATE> strmat(cmesh);
//         strmat.SetNumThreads(0); // ou outro número
//
//
//         analysis.SetStructuralMatrix(strmat);
//
//         analysis.StructMatrix()->EquationFilter().Reset();
//
//         int nact = analysis.StructMatrix()->EquationFilter().NActiveEquations();
//
//         TPZKrylovEigenSolver<STATE> solver;
//
//
//         int nev  = std::min(24, nact-1); // número de autovalores
//
//         solver.SetNEigenpairs(nev);
//
//         int kdim = std::min(nact, std::max(nev+10, 4*nev));
//
//         solver.SetKrylovDim(kdim);
//
//         solver.SetTolerance(1e-12);
//
//         solver.SetAsGeneralised(true);
//
//         analysis.SetSolver(solver);
//         analysis.Assemble();
//         analysis.Solve();
//
//         cout << "f" << endl;
//         TPZVec<CSTATE> vec=analysis.GetEigenvalues();
//
//         for(int i=0;i<vec.size();i++)cout<< vec[i].real() <<endl;
//
//
//         return 0;
// }
