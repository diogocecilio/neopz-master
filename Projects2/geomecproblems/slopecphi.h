#pragma once

#include "pzcmesh.h"

#include "pzgmesh.h"
#include "TPZYCMohrCoulombPV.h"

#include "TPZPlasticStepPV.h"
#include "pzelastoplasticanalysis.h"
#include <pzskylstrmatrix.h> //symmetric skyline matrix storage
#include <pzstepsolver.h> //for TPZStepSolver
#include "pzfstrmatrix.h"
#include "readgidmesh.h"

typedef TPZPlasticStepPV<TPZYCMohrCoulombPV, TPZElasticResponse> LEMC;
//typedef TPZPlasticStepPV<TPZYCMohrCoulombPV, TPZElasticResponse> LEVM;

typedef   TPZMatElastoPlastic2D <LEMC, TPZElastoPlasticMem > plasticmat;

class SlopeCphi
{
protected:



public:

TPZGeoMesh *CreateGeoMesh();
TPZCompMesh *CreateCompMesh(TPZGeoMesh *gmesh);
void PostProcessVariables(TPZStack<std::string> &scalNames, TPZStack<std::string> &vecNames);
void  CreatePostProcessingMesh (TPZPostProcAnalysis * PostProcess ,TPZCompMesh * cmesh);
void PostPlasticity(string vtkd,TPZCompMesh * cmesh);
TPZVec<REAL> findid(TPZCompMesh *cmesh,TPZVec<REAL> coord);
void run();
void ApplyGravityLoad(TPZManVector<REAL, 3> bodyforce);
};


void SlopeCphi::run(){


    TPZGeoMesh *gmesh =  CreateGeoMesh();

     TPZCompMesh*cmesh = CreateCompMesh(gmesh);
//
     TPZElastoPlasticAnalysis * analysis = new TPZElastoPlasticAnalysis(cmesh,cout);
//
//
    TPZManVector<REAL,3> gravityload(3,0);
    gravityload[1]=-20.;
    plasticmat * body= dynamic_cast<plasticmat *> ( cmesh->FindMaterial ( 1 ) );
    body->SetBodyForce ( gravityload );

     TPZSkylineStructMatrix<REAL> matskl ( cmesh );

     matskl.SetNumThreads ( 6 );
//
     TPZStepSolver<STATE> step;
//
     step.SetDirect ( ELDLt );
//
     analysis->SetStructuralMatrix(matskl);
//
     analysis->SetSolver ( step );
//
     REAL tolfs=0.001;
     REAL numiterfs=200;
     REAL numiterres =20;
     REAL tolres =1.e-5;
     REAL l =0.5;
     REAL lambda0=0.2;
     bool converge;
     REAL ndesirediters=8;
     REAL llimit=0.1;
REAL tol=1.e-3;
 	int numiter=10;
 	bool linesearch=true;
 	bool checkconv= false;
    bool convordiv;
    int iters;
     analysis->IterativeProcess(std::cout,tol,numiter,linesearch,checkconv);
  //   REAL fs = analysis->IterativeProcessHybridArcLength ( tolfs,numiterfs,tolres,numiterres,l,lambda0,converge,ndesirediters,llimit );
    // REAL fs = analysis->IterativeProcessLinearisedArcLength()
     //REAL fs = analysis->IterativeProcessArcLength ( tolfs,numiterfs,tolres,numiterres,l,lambda0,converge );
     string file  = "slopecphi.vtk";
     PostPlasticity(file,cmesh);


}

TPZGeoMesh *SlopeCphi::CreateGeoMesh() {

    string file = "/home/diogo/projects/pz-random/data/meshcphi.msh";

    TPZGeoMesh *gmesh  =  new TPZGeoMesh();

    gmesh->SetDimension ( 2 );

    readgidmesh read;

    std::vector<std::vector<int>> meshtopol;
    std::vector<std::vector<double>> meshcoords;

    read.ReadMesh2 ( meshtopol,meshcoords,file );


    cout << "a" << endl;
    int ncoords = meshcoords.size();
    gmesh->NodeVec().Resize ( ncoords );

    TPZVec<REAL> coord ( 2 );
    for ( int inode=0; inode<ncoords; inode++ )
    {
        coord[0] = meshcoords[inode][1];
        coord[1] = meshcoords[inode][2];
        gmesh->NodeVec() [inode] = TPZGeoNode ( inode, coord, *gmesh );
    }

    int sz = meshtopol.size();
    TPZVec <long> TopoQuad ( 4 );
    TPZVec <long> TopoTri ( 3 );
    TPZVec <long> TopoLine ( 2 );
    for ( int iel=0; iel<meshtopol.size(); iel++ )
    {
        if ( meshtopol[iel].size() ==4 )
        {
            TopoTri[0] =meshtopol[iel][1]-1;
            TopoTri[1] =meshtopol[iel][2]-1;
            TopoTri[2] =meshtopol[iel][3]-1;
            new TPZGeoElRefPattern< pzgeom::TPZGeoTriangle> ( iel, TopoTri, 1,*gmesh );

        }
        else if ( meshtopol[iel].size() ==5 )
        {

            TopoQuad[0] =meshtopol[iel][1]-1;
            TopoQuad[1] =meshtopol[iel][2]-1;
            TopoQuad[2] =meshtopol[iel][3]-1;
            TopoQuad[3] =meshtopol[iel][4]-1;
            new TPZGeoElRefPattern< pzgeom::TPZGeoQuad> ( iel, TopoQuad, 1,*gmesh );

        }
        else if ( meshtopol[iel].size() ==3 )
        {
            TopoLine[0] = meshtopol[iel][1]-1;
            TopoLine[1] = meshtopol[iel][2]-1;
            REAL x0 = meshcoords[TopoLine[0]][1];
            REAL y0 = meshcoords[TopoLine[0]][2];
            REAL xf = meshcoords[TopoLine[1]][1];
            REAL yf = meshcoords[TopoLine[1]][2];
            REAL tol=1.e-3;
            if ( ( fabs ( ( y0-0 ) ) <tol && fabs ( ( yf-0 ) ) <tol ) )
            {
                new TPZGeoElRefPattern< pzgeom::TPZGeoLinear> ( iel, TopoLine, -1, *gmesh );//bottom
            }
            else if ( ( fabs ( ( x0-75 ) ) <tol && fabs ( ( xf-75 ) ) <tol ) )
            {
                new TPZGeoElRefPattern< pzgeom::TPZGeoLinear> ( iel, TopoLine, -2, *gmesh );//rigth
            }
            else if ( ( fabs ( ( y0-30 ) ) <tol && fabs ( ( yf-30 ) ) <tol ) )
            {
                new TPZGeoElRefPattern< pzgeom::TPZGeoLinear> ( iel, TopoLine, -3, *gmesh );//toprigth
            }
            else if ( ( fabs ( ( y0-40 ) ) <tol && fabs ( ( yf-40 ) ) <tol ) )
            {
                new TPZGeoElRefPattern< pzgeom::TPZGeoLinear> ( iel, TopoLine, -4, *gmesh );//topleft
            }
            else if ( ( fabs ( ( x0-0 ) ) <tol && fabs ( ( xf-0 ) ) <tol ) )
            {
                new TPZGeoElRefPattern< pzgeom::TPZGeoLinear> ( iel, TopoLine, -5, *gmesh );//left
            }
            else if ( ( fabs ( ( xf-x0 ) ) >tol && fabs ( ( yf-y0 ) ) >tol ) )
            {
                new TPZGeoElRefPattern< pzgeom::TPZGeoLinear> ( iel, TopoLine, -6, *gmesh );//ramp
            }
            else
            {
                cout<< "bc element not found."<<endl;
                cout<< "x0 = " << x0 << " y0 = "<< y0 << endl;
                cout<< "xf = " << xf << " yf = "<< yf << endl;
               // DebugStop();
            }

        }

    }


    gmesh->BuildConnectivity();
    cout << "c" << endl;
    for ( int d = 0; d<3; d++ )
    {
        int nel = gmesh->NElements();
        TPZManVector<TPZGeoEl *> subels;
        for ( int iel = 0; iel<nel; iel++ )
        {
            TPZGeoEl *gel = gmesh->ElementVec() [iel];
            gel->Divide ( subels );
        }
    }
    std::ofstream print("gmeshslopecphi.txt");
	gmesh->Print(print);
    std::ofstream files ( "teste-mesh-slopecphi.vtk" );
    TPZVTKGeoMesh::PrintGMeshVTK ( gmesh,files,false );
    cout << "d" << endl;
    return gmesh;

}
//     REAL co[14][2] = {
//         {75,0}, {75,30},{37.5,0},{40,25},
//     {40,30},{35,25},{35,30},{30,25},
//     {30,30},{35,40},{30,40},{0,0},{0,30},{0,40}};
//
// long indices[18][4] = {{14 ,13, 9, 11},{12, 8 ,9 ,13},{10 ,11, 9 ,7},{7 ,6, 4, 5},
// {9 ,8, 6, 7},{2 ,5 ,4 ,1},{3, 6, 8, 12},{1 ,4, 6, 3},{3, 1},{12 ,3},{2 ,5},{1 ,2},{5 ,7},{10,11},{7 ,10},{13, 12},{11 ,14},{14 ,13}
// };
void SlopeCphi::ApplyGravityLoad(TPZManVector<REAL, 3> bodyforce)
{
   // plasticmat * body= dynamic_cast<plasticmat *> ( fCompMesh->FindMaterial ( 1 ) );
   // body->SetBodyForce ( bodyforce );

}
TPZCompMesh *SlopeCphi::CreateCompMesh(TPZGeoMesh *gmesh) {


    unsigned int dim  = 2;
    const std::string name ( "ElastoPlastic COMP MESH Slope Problem " );

    TPZCompMesh * cmesh =  new TPZCompMesh ( gmesh );
    cmesh->SetName ( name );
    cmesh->SetDefaultOrder (2 );
    cmesh->SetDimModel ( dim );


    //Pu = (2+pi)c = 218.183

      // Mohr Coulomb data
    //REAL mc_cohesion    = 490.;//kPa
    REAL mc_cohesion    = 50.;//N/cm^2
    REAL mc_phi         = 20.*M_PI/180;
    REAL mc_psi         = mc_phi;

    /// ElastoPlastic Material using Mohr Coulomb
    // Elastic predictor
    TPZElasticResponse ER;
    REAL nu = 0.49;
    //REAL E = 10000000.;////kPa
    REAL E = 20000.;////N/cm^2

    LEMC plasticstep;
    //TPZPlasticStepPV<TPZYCMohrCoulombPV, TPZElasticResponse> LEMC;
    //TPZPlasticStepVoigt<TPZMohrCoulombVoigt, TPZElasticResponse> LEMC;
    ER.SetEngineeringData(E,nu);
    //ER.SetUp ( E, nu );
    plasticstep.fER =ER;
    // LEMC.SetElasticResponse( ER );
    plasticstep.fYC.SetUp ( mc_phi, mc_psi, mc_cohesion, ER );

    int PlaneStrain = 1;

    plasticmat * material = new plasticmat ( 1,PlaneStrain );

    REAL factor;
    TPZManVector<REAL, 3> bodyforce ( 3,0. );
    factor=1.;
    bodyforce[1]=-20.;

   // material->SetPlasticity ( plasticstep );

    material->SetId ( 1 );

   // material->SetWhichLoadVector ( 0 ); //option to compute the total internal force vecor fi=(Bt sigma+ N (b+gradu))

   // material->SetLoadFactor ( factor );

    material->SetBodyForce ( bodyforce );

    cmesh->InsertMaterialObject(material);

    TPZFMatrix<STATE>  val1 ( 2,2,0. );

    TPZManVector<STATE,2> val2 ( 2,0. );

    val2[1]=1.;
    val2[0]=1.;
    auto * bcleft = material->CreateBC(material,-1,3,val1,val2);//bottom

    val2[0]=1.;
    val2[1]=0.;
    auto * bcbottom = material->CreateBC(material,-2,3,val1,val2);//rigth

    val2[0]=1.;
    val2[1]=0;
    auto *  bcrigth = material->CreateBC(material,-3,3,val1,val2);//left


	cmesh->InsertMaterialObject(bcleft);

    cmesh->InsertMaterialObject(bcbottom);

    cmesh->InsertMaterialObject(bcrigth);



	cmesh->SetAllCreateFunctionsContinuousWithMem();

    cmesh->AutoBuild();

    std::ofstream print("cmeshslopecphi.txt");
	cmesh->Print(print);

    return cmesh;
}

void SlopeCphi::PostProcessVariables ( TPZStack<std::string> &scalNames, TPZStack<std::string> &vecNames )
{




    scalNames.Push ( "PlStrainSqJ2" );
    scalNames.Push ( "PlStrainSqJ2El" );
    scalNames.Push ( "PlAlpha" );

    scalNames.Push ( "DisplacementX" );
    scalNames.Push ( "DisplacementY" );
    scalNames.Push ( "DisplacementZ" );
    vecNames.Push ( "DisplacementTotal" );
    vecNames.Push ( "PrincipalStress" );


}

void SlopeCphi::PostPlasticity(string vtkd,TPZCompMesh * cmesh)
{
    TPZPostProcAnalysis * postprocdeter = new TPZPostProcAnalysis();
    CreatePostProcessingMesh ( postprocdeter ,cmesh);

    TPZVec<int> PostProcMatIds ( 1,1 );

    TPZStack<std::string> PostProcVars, scalNames, vecNames;

    PostProcessVariables ( scalNames, vecNames );

    //string vtkd = "postprocessdeter.vtk";
    postprocdeter->DefineGraphMesh ( 2,scalNames,vecNames,vtkd );

    postprocdeter->PostProcess ( 0 );

    delete postprocdeter;
}
void  SlopeCphi::CreatePostProcessingMesh (TPZPostProcAnalysis * PostProcess ,TPZCompMesh * cmesh)
{
    if ( PostProcess->ReferenceCompMesh() != cmesh )
    {

        PostProcess->SetCompMesh ( cmesh );

        TPZVec<int> PostProcMatIds ( 1,1 );
        TPZStack<std::string> PostProcVars, scalNames, vecNames;
        PostProcessVariables ( scalNames, vecNames );

        for ( int i=0; i<scalNames.size(); i++ )
        {
            PostProcVars.Push ( scalNames[i] );
        }
        for ( int i=0; i<vecNames.size(); i++ )
        {
            PostProcVars.Push ( vecNames[i] );
        }
        //
        PostProcess->SetPostProcessVariables ( PostProcMatIds, PostProcVars );
        TPZFStructMatrix<REAL> structmatrix ( PostProcess->Mesh() );
        structmatrix.SetNumThreads ( 0 );
        PostProcess->SetStructuralMatrix ( structmatrix );
    }
    //
    //Chamar com o analysis e nao com o postanalysis pois tem o acumulo de sols
    PostProcess->TransferSolution();

}

TPZVec<REAL> SlopeCphi::findid(TPZCompMesh *cmesh,TPZVec<REAL> coord) {

    TPZGeoMesh * gmesh =  cmesh->Reference();
    int dim = gmesh->Dimension()+1;
    TPZVec<REAL> xd(dim,0.);
    TPZVec<REAL> mpt(dim,0.);
    xd[0] =coord[0];
    xd[1] = coord[1];
    int id;
    int nels = gmesh->NElements();
    for(int iel=0; iel<nels; iel++) {


        TPZVec<REAL> qsi(2,0.);
        TPZGeoEl * gel = gmesh->ElementVec()[iel];

        TPZCompEl * cel = gel->Reference();



        if(gel->MaterialId()<0)continue;
        bool check = gel->ComputeXInverse(xd,qsi,1.e-5);

        if(check==true)
        {

            int nnodes=gel->NNodes();
            for(int inode=0; inode<nnodes; inode++)
            {
                TPZVec<REAL> co(dim);
                TPZGeoNode* node = gel->NodePtr(inode);
                node->GetCoordinates(co);
                id = node->Id();

                if(fabs(co[0]-xd[0])<1.e-3 && fabs(co[1]-xd[1])<1.e-3)
                {
                    TPZInterpolationSpace *intel = dynamic_cast<TPZInterpolationSpace *> ( cel );

//                     TPZMaterialData data;
//
//                     data.fNeedsSol = true;
//
//                     intel->InitMaterialData ( data );
//
//                     intel->ComputeRequiredData ( data, qsi );
//
//  					//cout << " data.sol  = "<<data.sol     << endl;
//                     return data.sol[0];
                }
            }


        }
    }

}
