#include <iostream>
#include <fstream>
#include <thread>
#include <vector>
#include "TPZFileStream.h"
#include <TPZBFileStream.h>
#include "SlopeAnalysis.h"
#include <fstream>
#include <iostream>
#include <fstream>
#include <thread>
#include <vector>
#include <mutex>
#include <iostream>
#include <fstream>
#include <unistd.h>
#include <sys/wait.h>
#include <memory>
#include "SlopeAnalysis.h"
#include "TPZEigenSolver.h"
#include "Plasticity/TPZElasticResponse.h"
#include "TPZKrylovEigenSolver.h"
#include "TPZLapackEigenSolver.h" // ou outro solver concreto
typedef TPZPlasticStepPV<TPZYCMohrCoulombPV, TPZElasticResponse> plasticmorh;
typedef TPZMatElastoPlastic2D <TPZPlasticStepPV<TPZYCMohrCoulombPV, TPZElasticResponse>, TPZElastoPlasticMem> plasticmat;

TPZGeoMesh * TriGMesh ( int ref );

TPZCompMesh * CreateCMesh ( TPZGeoMesh *gmesh, int pOrder, plasticmat * mat);

plasticmat * CreateMaterial(REAL young, REAL poisson, REAL coes,REAL atrito,TPZManVector<REAL, 3> bodyforce,int planestrain=1,int matid=1);



int main()
{

        std::cout << "\ncriando matrizes:\n";
        TPZFMatrix<REAL> fullmatA(2,2);
        fullmatA.PutVal(0,0,12.);fullmatA.PutVal(0,1,-4.);
        fullmatA.PutVal(1,0,-6.);fullmatA.PutVal(1,1,8.);
        fullmatA.Print(cout);
        TPZFMatrix<REAL> fullmatB(2,2);
        fullmatB.PutVal(0,0,7.);fullmatB.PutVal(0,1,-1.);
        fullmatB.PutVal(1,0,-9.);fullmatB.PutVal(1,1,5.);
        fullmatB.Print(cout);

        // Criar ponteiro para a matriz
        TPZAutoPointer<TPZFMatrix<REAL>> A = new TPZFMatrix<REAL>(fullmatA);



        // std::cout << "\nadicionado matrizes\n";
        // eigensolver->SetMatrixA(&fullmatA);
        // eigensolver->SetMatrixB(&fullmatB);
        // TPZVec<std::complex<REAL>> val;
        // TPZFMatrix<std::complex<REAL>> vec;
        // std::cout << "\nresolvendo\n";
        // eigensolver->SolveGeneralisedEigenProblem(val,vec);
        //
        // std::cout << "\nAutovalores:\n";
        // for (auto &v : val) std::cout << v << "\n";
        //
        // std::cout << "\nAutovetores:\n";
        // vec.Print(std::cout);

        return 0;



        int ref =1;
        //TPZGeoMesh * gmesh =  TriGMesh (ref);
        auto gmesh =  TriGMesh (ref);

        REAL young=20000.;
        REAL poisson=0.2;
        REAL coes=10.;
        REAL atrito=30*M_PI/180.;
        TPZManVector<REAL, 3> bodyforce ( 3,0. );
        bodyforce[1]=-20.;
        plasticmat * mat = CreateMaterial(young,poisson,coes,atrito,bodyforce);

        int pOrder=2;
        TPZCompMesh * cmesh = CreateCMesh (gmesh, pOrder, mat);

        SlopeAnalysis * plasticanal= new SlopeAnalysis(gmesh,cmesh);

        bool issrm=true;
        plasticanal->SolveDeterministic(issrm,coes,atrito);

        string vtkd = "post.vtk";
        plasticanal->PostPlasticity (vtkd );


        return 0;
}

TPZGeoMesh * TriGMesh ( int ref )
{

        TPZGeoMesh *gmesh  =  new TPZGeoMesh();


        gmesh->SetDimension ( 2 );

        TPZVec<REAL> coord ( 2 );

        vector<vector<double>> co= {
                /*0*/{0,0},/*1*/{10,0},/*2*/{20,0},/*3*/{30,0},/*4*/{40,0},/*5*/{50,0},/*6*/{60,0},/*7*/{70,0},
                /*8*/{0,10},/*9*/{10,10},/*10*/{20,10},/*11*/{30,10},/*12*/{40,10},/*13*/{50,10},/*14*/{60,10},/*15*/{70,10},
                /*16*/{0,20},/*17*/{10,20},/*18*/{20,20},/*19*/{30,20},/*20*/{40,20},/*21*/{50,20},/*22*/{60,20},/*23*/{70,20},
                /*24*/{0,30},/*25*/{10,30},/*26*/{20,30},/*27*/{30,30},/*28*/{40,30},/*29*/{50,30},/*30*/{60,30},/*31*/{70,30},
                /*32*/{0,40},/*33*/{10,40},/*34*/{20,40},/*35*/{30,40}
        };
        vector<vector<int>> topol = {

                /*0*/{0,1,8},/*1*/{1,9,8},/*2*/{1,2,9},/*3*/{2,10,9},/*4*/{2,3,10},/*5*/{3,11,10},/*6*/{3,4,11},
                /*7*/{4,12,11},/*8*/{4,5,12}/*9*/,{5,13,12},/*10*/{5,6,13},/*11*/{6,14,13},/*12*/{6,7,14},/*13*/{7,15,14},

                /*14*/{8,9,16},/*15*/{9,17,16},/*16*/{9,10,17},/*17*/{10,18,17},/*18*/{10,11,18},/*19*/{11,19,18},/*20*/{11,12,19},
                /*21*/{12,20,19},/*22*/{12,13,20}/*23*/,{13,21,20},/*24*/{13,14,21},/*25*/{14,22,21},/*26*/{14,15,22},/*27*/{15,23,22},

                /*28*/{16,17,24},/*29*/{17,25,24},/*30*/{17,18,25},/*31*/{18,26,25},/*32*/{18,19,26},/*33*/{19,27,26},/*34*/{19,20,27},
                /*35*/{20,28,27},/*36*/{20,21,28}/*37*/,{21,29,28},/*38*/{21,22,29},/*39*/{22,30,29},/*40*/{22,23,30},/*41*/{23,31,30},

                /*42*/{24,25,32},/*43*/{25,33,32},/*44*/{25,26,33},/*45*/{26,34,33},/*46*/{26,27,34},/*47*/{27,35,34},/*48*/{27,28,35},

                {0,1},{1,2},{2,3},{3,4},{4,5},{5,6},{6,7},/*-1 bottom*/

                {7,15},{15,23},{23,31},/*-2 right*/

                {31,30},{30,29},{29,28},/*-3 top right*/

                {35,34},{34,33},{33,32},/*-4 top left*/

                {32,24},{24,16},{16,8},{8,0},/*-5 left*/

                {28,35}/*-6 ramp*/


        };

        gmesh->NodeVec().Resize ( co.size() );

        for ( int inode=0; inode<co.size(); inode++ ) {
                coord[0] = co[inode][0];
                coord[1] = co[inode][1];
                gmesh->NodeVec() [inode] = TPZGeoNode ( inode, coord, *gmesh );
        }
        TPZVec <long> topotri ( 3 );
        TPZVec <long> TopoLine ( 2 );
        for ( int iel=0; iel<topol.size(); iel++ ) {
                if ( topol[iel].size() ==3 ) {
                        topotri[0] = topol[iel][0];
                        topotri[1] = topol[iel][1];
                        topotri[2] = topol[iel][2];
                        new TPZGeoElRefPattern< pzgeom::TPZGeoTriangle> ( iel, topotri, 1,*gmesh );
                } else if ( topol[iel].size() ==2 ) {

                        TopoLine[0] = topol[iel][0];
                        TopoLine[1] = topol[iel][1];
                        REAL x0 = co[TopoLine[0]][0];
                        REAL y0 = co[TopoLine[0]][1];
                        REAL xf = co[TopoLine[1]][0];
                        REAL yf = co[TopoLine[1]][1];
                        REAL tol=1.e-3;
                        REAL L=70;
                        REAL h1=30;
                        REAL h2=10;
                        if ( ( fabs ( ( y0-0 ) ) <tol && fabs ( ( yf-0 ) ) <tol ) ) {
                                //bottom
                                new TPZGeoElRefPattern< pzgeom::TPZGeoLinear> ( iel, TopoLine, -1, *gmesh );
                        } else if ( ( fabs ( ( x0-L ) ) <tol && fabs ( ( xf-L ) ) <tol ) ) {
                                //rigth
                                new TPZGeoElRefPattern< pzgeom::TPZGeoLinear> ( iel, TopoLine, -2, *gmesh );
                        } else if ( ( fabs ( ( y0-h1 ) ) <tol && fabs ( ( yf-h1 ) ) <tol ) ) {
                                //toprigth
                                new TPZGeoElRefPattern< pzgeom::TPZGeoLinear> ( iel, TopoLine, -3, *gmesh );
                        } else if ( ( fabs ( ( y0- ( h1+h2 ) ) ) <tol && fabs ( ( yf- ( h1+h2 ) ) ) <tol ) ) {
                                //topleft
                                new TPZGeoElRefPattern< pzgeom::TPZGeoLinear> ( iel, TopoLine, -4, *gmesh );
                        } else if ( ( fabs ( ( x0-0 ) ) <tol && fabs ( ( xf-0 ) ) <tol ) ) {
                                //left
                                new TPZGeoElRefPattern< pzgeom::TPZGeoLinear> ( iel, TopoLine, -5, *gmesh );
                        } else if ( ( fabs ( ( xf-x0 ) ) >tol && fabs ( ( yf-y0 ) ) >tol ) ) {
                                //ramp
                                new TPZGeoElRefPattern< pzgeom::TPZGeoLinear> ( iel, TopoLine, -6, *gmesh );
                        } else {
                                cout<< "bc element not found."<<endl;
                                cout<< "x0 = " << x0 << " y0 = "<< y0 << endl;
                                cout<< "xf = " << xf << " yf = "<< yf << endl;
                                DebugStop();
                        }

                }
        }

        gmesh->BuildConnectivity();
        for ( int d = 0; d<ref; d++ ) {
                int nel = gmesh->NElements();
                TPZManVector<TPZGeoEl *> subels;
                for ( int iel = 0; iel<nel; iel++ ) {
                        TPZGeoEl *gel = gmesh->ElementVec() [iel];
                        gel->Divide ( subels );
                }
        }

        string meshref = "gmeshtri.vtk";
        std::ofstream files ( meshref );
        TPZVTKGeoMesh::PrintGMeshVTK ( gmesh,files,true );
        return gmesh;
}

plasticmat * CreateMaterial(REAL young, REAL poisson, REAL coes,REAL atrito,TPZManVector<REAL, 3> bodyforce,int planestrain,int matid)
{

        TPZElasticResponse  elasticresponse;

        elasticresponse.SetEngineeringData (young,poisson );
        // Mohr Coulomb data
        REAL mc_cohesion    = coes;                         //kPa
        REAL mc_phi         = atrito;
        REAL mc_psi         = atrito;

        //elasticresponse.Print(std::cout);

        plasticmorh mohrcoulombplasticstep;

        mohrcoulombplasticstep.fYC.SetUp ( mc_phi, mc_psi, mc_cohesion, elasticresponse );

        mohrcoulombplasticstep.fER = elasticresponse;

        mohrcoulombplasticstep.SetStrengthReductionFactor(1.);

        plasticmat * material = new plasticmat ( matid,planestrain );

        material->SetPlasticityModel ( mohrcoulombplasticstep );

        material->SetId ( matid );

        REAL factor;
        factor=1.;

        material->SetId ( 1 );

        DebugStop();
        //REIMPLEMENTAR
        // material->SetWhichLoadVector ( 0 );//option to compute the total internal force vecor fi=(Bt sigma+ N (b+gradu))
        //
        // material->SetLoadFactor ( factor );
        //
        // material->SetBodyForce ( bodyforce );

        return material;
}

TPZCompMesh * CreateCMesh ( TPZGeoMesh *gmesh, int pOrder, plasticmat * mat)
{
        // Creating computational mesh:
        TPZCompMesh * cmesh = new TPZCompMesh ( gmesh );

        //cmesh->Print(std::cout);
        cmesh->SetDefaultOrder ( pOrder );

        int dim = 2 ;

        cmesh->SetDimModel ( dim );

        cmesh->InsertMaterialObject ( mat );

        //cmesh->Print(std::cout);

        // boundary condition
        TPZFMatrix<STATE>  val1 ( 2,2,0. );

        TPZManVector<STATE,2> val2 ( 2,0. );

        int directionaldirichlet = 3 ;

        val2[0]=1;
        val2[1]=1;
        auto * BCond0 = mat->CreateBC ( mat, -1, directionaldirichlet, val1, val2 );

        val2[0]=1;
        val2[1]=0;
        auto * BCond1 = mat->CreateBC ( mat, -2, directionaldirichlet, val1, val2 );

        val2[0]=1;
        val2[1]=0;
        auto * BCond2 = mat->CreateBC ( mat, -5, directionaldirichlet, val1, val2 );

        cmesh->InsertMaterialObject ( BCond0 );
        cmesh->InsertMaterialObject ( BCond1 );
        cmesh->InsertMaterialObject ( BCond2 );

        cmesh->SetAllCreateFunctionsContinuousWithMem();
        //Creating computational elements that manage the space of the mesh:
        cmesh->AutoBuild();
        cmesh->AdjustBoundaryElements();
        cmesh->CleanUpUnconnectedNodes();

        return cmesh;
}

// #include <iostream>
// #include <fstream>
// #include <vector>
// #include <memory>
// #include <cmath>
// #include "TPZFileStream.h"
// #include "TPZBFileStream.h"
// #include "SlopeAnalysis.h"
//
// // Namespace para evitar poluição global
// namespace GeomecProbDeter {
//
//         typedef TPZPlasticStepPV<TPZYCMohrCoulombPV, TPZElasticResponse> PlasticMohr;
//         typedef TPZMatElastoPlastic2D<PlasticMohr, TPZElastoPlasticMem> PlasticMat;
//
//         struct MaterialConfig {
//                 REAL young;
//                 REAL poisson;
//                 REAL coes;
//                 REAL atrito;
//                 TPZManVector<REAL, 3> bodyforce;
//                 int planestrain = 1;
//                 int matid = 1;
//         };
//
//         /**
//          * Cria a malha geométrica triangular.
//          * @param ref Quantidade de refinamentos
//          * @return ponteiro inteligente para TPZGeoMesh
//          */
//         std::unique_ptr<TPZGeoMesh> TriGMesh(int ref);
//
//         /**
//          * Cria o material elastoplástico
//          * @param config Configuração do material
//          * @return ponteiro inteligente para PlasticMat
//          */
//         std::unique_ptr<PlasticMat> CreateMaterial(const MaterialConfig& config);
//
//         /**
//          * Cria a malha computacional
//          * @param gmesh Ponteiro para a malha geométrica
//          * @param pOrder Ordem do polinômio
//          * @param mat Ponteiro para o material
//          * @return ponteiro inteligente para TPZCompMesh
//          */
//         std::unique_ptr<TPZCompMesh> CreateCMesh(TPZGeoMesh* gmesh, int pOrder, PlasticMat* mat);
//
// } // namespace GeomecProbDeter
//
// int main() {
//         using namespace GeomecProbDeter;
//
//         int ref = 1;
//         auto gmesh = TriGMesh(ref);
//
//         MaterialConfig config;
//         config.young = 20000.;
//         config.poisson = 0.2;
//         config.coes = 10.;
//         config.atrito = 30 * M_PI / 180.;
//         config.bodyforce = TPZManVector<REAL, 3>(3, 0.);
//         config.bodyforce[1] = -20.;
//
//         auto mat = CreateMaterial(config);
//
//         int pOrder = 2;
//         auto cmesh = CreateCMesh(gmesh.get(), pOrder, mat.get());
//
//         auto plasticanal = std::make_unique<SlopeAnalysis>(gmesh.get(), cmesh.get());
//
//         bool issrm = true;
//         plasticanal->SolveDeterministic(issrm, config.coes, config.atrito);
//
//         string vtkd = "post.vtk";
//         plasticanal->PostPlasticity (vtkd );
//
//         cout<< "FIM"<<endl;
//
//         return 0;
// }
//
// /**
//  * Cria a malha geométrica triangular.
//  */
// std::unique_ptr<TPZGeoMesh> GeomecProbDeter::TriGMesh(int ref)
// {
//         auto gmesh = std::make_unique<TPZGeoMesh>();
//         gmesh->SetDimension(2);
//         TPZVec<REAL> coord(2);
//
//         std::vector<std::vector<double>> co = {
//                 /*0*/{0,0},/*1*/{10,0},/*2*/{20,0},/*3*/{30,0},/*4*/{40,0},/*5*/{50,0},/*6*/{60,0},/*7*/{70,0},
//                 /*8*/{0,10},/*9*/{10,10},/*10*/{20,10},/*11*/{30,10},/*12*/{40,10},/*13*/{50,10},/*14*/{60,10},/*15*/{70,10},
//                 /*16*/{0,20},/*17*/{10,20},/*18*/{20,20},/*19*/{30,20},/*20*/{40,20},/*21*/{50,20},/*22*/{60,20},/*23*/{70,20},
//                 /*24*/{0,30},/*25*/{10,30},/*26*/{20,30},/*27*/{30,30},/*28*/{40,30},/*29*/{50,30},/*30*/{60,30},/*31*/{70,30},
//                 /*32*/{0,40},/*33*/{10,40},/*34*/{20,40},/*35*/{30,40}
//         };
//         std::vector<std::vector<int>> topol = {
//                 /*0*/{0,1,8},/*1*/{1,9,8},/*2*/{1,2,9},/*3*/{2,10,9},/*4*/{2,3,10},/*5*/{3,11,10},/*6*/{3,4,11},
//                 /*7*/{4,12,11},/*8*/{4,5,12},/*9*/{5,13,12},/*10*/{5,6,13},/*11*/{6,14,13},/*12*/{6,7,14},/*13*/{7,15,14},
//
//                 /*14*/{8,9,16},/*15*/{9,17,16},/*16*/{9,10,17},/*17*/{10,18,17},/*18*/{10,11,18},/*19*/{11,19,18},/*20*/{11,12,19},
//                 /*21*/{12,20,19},/*22*/{12,13,20},/*23*/{13,21,20},/*24*/{13,14,21},/*25*/{14,22,21},/*26*/{14,15,22},/*27*/{15,23,22},
//
//                 /*28*/{16,17,24},/*29*/{17,25,24},/*30*/{17,18,25},/*31*/{18,26,25},/*32*/{18,19,26},/*33*/{19,27,26},/*34*/{19,20,27},
//                 /*35*/{20,28,27},/*36*/{20,21,28},/*37*/{21,29,28},/*38*/{21,22,29},/*39*/{22,30,29},/*40*/{22,23,30},/*41*/{23,31,30},
//
//                 /*42*/{24,25,32},/*43*/{25,33,32},/*44*/{25,26,33},/*45*/{26,34,33},/*46*/{26,27,34},/*47*/{27,35,34},/*48*/{27,28,35},
//
//                 {0,1},{1,2},{2,3},{3,4},{4,5},{5,6},{6,7},/*-1 bottom*/
//
//                 {7,15},{15,23},{23,31},/*-2 right*/
//
//                 {31,30},{30,29},{29,28},/*-3 top right*/
//
//                 {35,34},{34,33},{33,32},/*-4 top left*/
//
//                 {32,24},{24,16},{16,8},{8,0},/*-5 left*/
//
//                 {28,35}/*-6 ramp*/
//         };
//
//         gmesh->NodeVec().Resize(co.size());
//         for (size_t inode = 0; inode < co.size(); inode++) {
//                 coord[0] = co[inode][0];
//                 coord[1] = co[inode][1];
//                 gmesh->NodeVec()[inode] = TPZGeoNode(inode, coord, *gmesh);
//         }
//
//         TPZVec<long> topotri(3);
//         TPZVec<long> TopoLine(2);
//         for (size_t iel = 0; iel < topol.size(); iel++) {
//                 if (topol[iel].size() == 3) {
//                         topotri[0] = topol[iel][0];
//                         topotri[1] = topol[iel][1];
//                         topotri[2] = topol[iel][2];
//                         new TPZGeoElRefPattern<pzgeom::TPZGeoTriangle>(iel, topotri, 1, *gmesh);
//                 } else if (topol[iel].size() == 2) {
//                         TopoLine[0] = topol[iel][0];
//                         TopoLine[1] = topol[iel][1];
//                         REAL x0 = co[TopoLine[0]][0];
//                         REAL y0 = co[TopoLine[0]][1];
//                         REAL xf = co[TopoLine[1]][0];
//                         REAL yf = co[TopoLine[1]][1];
//                         REAL tol = 1.e-3;
//                         REAL L = 70;
//                         REAL h1 = 30;
//                         REAL h2 = 10;
//
//                         if (std::fabs(y0 - 0) < tol && std::fabs(yf - 0) < tol) {
//                                 new TPZGeoElRefPattern<pzgeom::TPZGeoLinear>(iel, TopoLine, -1, *gmesh);
//                         } else if (std::fabs(x0 - L) < tol && std::fabs(xf - L) < tol) {
//                                 new TPZGeoElRefPattern<pzgeom::TPZGeoLinear>(iel, TopoLine, -2, *gmesh);
//                         } else if (std::fabs(y0 - h1) < tol && std::fabs(yf - h1) < tol) {
//                                 new TPZGeoElRefPattern<pzgeom::TPZGeoLinear>(iel, TopoLine, -3, *gmesh);
//                         } else if (std::fabs(y0 - (h1 + h2)) < tol && std::fabs(yf - (h1 + h2)) < tol) {
//                                 new TPZGeoElRefPattern<pzgeom::TPZGeoLinear>(iel, TopoLine, -4, *gmesh);
//                         } else if (std::fabs(x0 - 0) < tol && std::fabs(xf - 0) < tol) {
//                                 new TPZGeoElRefPattern<pzgeom::TPZGeoLinear>(iel, TopoLine, -5, *gmesh);
//                         } else if (std::fabs(xf - x0) > tol && std::fabs(yf - y0) > tol) {
//                                 new TPZGeoElRefPattern<pzgeom::TPZGeoLinear>(iel, TopoLine, -6, *gmesh);
//                         } else {
//                                 std::cerr << "bc element not found." << std::endl;
//                                 std::cerr << "x0 = " << x0 << " y0 = " << y0 << std::endl;
//                                 std::cerr << "xf = " << xf << " yf = " << yf << std::endl;
//                                 throw std::runtime_error("Boundary condition element not found.");
//                         }
//                 }
//         }
//
//         gmesh->BuildConnectivity();
//         for (int d = 0; d < ref; d++) {
//                 int nel = gmesh->NElements();
//                 TPZManVector<TPZGeoEl *> subels;
//                 for (int iel = 0; iel < nel; iel++) {
//                         TPZGeoEl *gel = gmesh->ElementVec()[iel];
//                         gel->Divide(subels);
//                 }
//         }
//
//         std::ofstream files("gmeshtri.vtk");
//         TPZVTKGeoMesh::PrintGMeshVTK(gmesh.get(), files, true);
//         return gmesh;
// }
//
// /**
//  * Cria o material elastoplástico
//  */
// std::unique_ptr<GeomecProbDeter::PlasticMat> GeomecProbDeter::CreateMaterial(const MaterialConfig& config)
// {
//         if (config.young <= 0 || config.poisson < 0 || config.coes < 0) {
//                 throw std::invalid_argument("Parâmetros físicos inválidos.");
//         }
//
//         TPZElasticResponse elasticresponse;
//         elasticresponse.SetEngineeringData(config.young, config.poisson);
//
//         PlasticMohr mohrcoulombplasticstep;
//         mohrcoulombplasticstep.fYC.SetUp(config.atrito, config.atrito, config.coes, elasticresponse);
//         mohrcoulombplasticstep.fER = elasticresponse;
//         mohrcoulombplasticstep.SetStrengthReductionFactor(1.);
//
//         auto material = std::make_unique<PlasticMat>(config.matid, config.planestrain);
//         material->SetPlasticityModel(mohrcoulombplasticstep);
//         material->SetId(config.matid);
//         material->SetWhichLoadVector(0);
//         material->SetLoadFactor(1.0);
//         material->SetBodyForce(config.bodyforce);
//
//         return material;
// }
//
// /**
//  * Cria a malha computacional
//  */
// std::unique_ptr<TPZCompMesh> GeomecProbDeter::CreateCMesh(TPZGeoMesh* gmesh, int pOrder, PlasticMat* mat)
// {
//         auto cmesh = std::make_unique<TPZCompMesh>(gmesh);
//         cmesh->SetDefaultOrder(pOrder);
//         cmesh->SetDimModel(2);
//         cmesh->InsertMaterialObject(mat);
//
//         TPZFMatrix<STATE> val1(2,2,0.);
//         TPZManVector<STATE,2> val2(2,0.);
//         int directionaldirichlet = 3;
//
//         val2[0]=1; val2[1]=1;
//         auto* BCond0 = mat->CreateBC(mat, -1, directionaldirichlet, val1, val2);
//
//         val2[0]=1; val2[1]=0;
//         auto* BCond1 = mat->CreateBC(mat, -2, directionaldirichlet, val1, val2);
//
//         val2[0]=1; val2[1]=0;
//         auto* BCond2 = mat->CreateBC(mat, -5, directionaldirichlet, val1, val2);
//
//         cmesh->InsertMaterialObject(BCond0);
//         cmesh->InsertMaterialObject(BCond1);
//         cmesh->InsertMaterialObject(BCond2);
//
//         cmesh->SetAllCreateFunctionsContinuousWithMem();
//         cmesh->AutoBuild();
//         cmesh->AdjustBoundaryElements();
//         cmesh->CleanUpUnconnectedNodes();
//
//         return cmesh;
// }
