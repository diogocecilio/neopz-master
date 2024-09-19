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

class Footing
{
protected:



public:

TPZGeoMesh *CreateGeoMesh();
TPZCompMesh *CreateCompMesh(TPZGeoMesh *gmesh);
void PostProcessVariables(TPZStack<std::string> &scalNames, TPZStack<std::string> &vecNames);
void  CreatePostProcessingMesh (TPZPostProcAnalysis * PostProcess ,TPZCompMesh * cmesh);
void PostPlasticity(string vtkd,TPZCompMesh * cmesh);
void findid(TPZCompMesh *cmesh,TPZVec<REAL> coord);
void runfoot();

};


void Footing::runfoot(){


    TPZGeoMesh *gmesh =  CreateGeoMesh();

     TPZCompMesh*cmesh = CreateCompMesh(gmesh);
//
     TPZElastoPlasticAnalysis * analysis = new TPZElastoPlasticAnalysis(cmesh,cout);
//
//
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
     REAL l =0.04;
     REAL lambda0=2;
     REAL ndesirediters=3;
     REAL llimit=0.04;
     bool converge;
REAL tol=1.e-3;
 	int numiter=10;
 	bool linesearch=true;
 	bool checkconv= false;
    bool convordiv;
    int iters;
     analysis->IterativeProcess(std::cout,tol,numiter,linesearch,checkconv);
    // REAL fs = analysis->IterativeProcessLinearisedArcLength()
     //REAL fs = analysis->IterativeProcessArcLength ( tolfs,numiterfs,tolres,numiterres,l,lambda0,converge );
     string file  = "footing.vtk";
     //PostPlasticity(file,cmesh);
 analysis->AcceptSolution();

     TPZLinearAnalysis * anal = new TPZLinearAnalysis(cmesh);
     //analysis->Solution().Print("sol");
    //TPZPostProcAnalysis * postprocdeter = new TPZPostProcAnalysis();
    //CreatePostProcessingMesh ( postprocdeter);
   // TPZNonLinearAnalysis * postprocdeter = new TPZNonLinearAnalysis(fCompMesh,cout);

    TPZVec<int> PostProcMatIds ( 1,1 );

    TPZStack<std::string> PostProcVars, scalNames, vecNames;

    PostProcessVariables ( scalNames, vecNames );

    //string vtkd = "postprocessdeter.vtk";
    anal->DefineGraphMesh ( 2,scalNames,vecNames,file );

    anal->PostProcess ( 0 );

}

TPZGeoMesh *Footing::CreateGeoMesh() {

    string file = "/home/diogo/projects/pz-random/data/mesh-foot.msh";

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
            else if ( ( fabs ( ( x0-500 ) ) <tol && fabs ( ( xf-500 ) ) <tol ) )
            {
                new TPZGeoElRefPattern< pzgeom::TPZGeoLinear> ( iel, TopoLine, -2, *gmesh );//rigth
            }
            else if ( ( fabs ( ( x0-0 ) ) <tol && fabs ( ( xf-0 ) ) <tol ) )
            {
                new TPZGeoElRefPattern< pzgeom::TPZGeoLinear> ( iel, TopoLine, -3, *gmesh );//left
            }
            else if ( (xf <=5   && fabs ( ( yf-500 ) ) <tol)  )
            {
                new TPZGeoElRefPattern< pzgeom::TPZGeoLinear> ( iel, TopoLine, -4, *gmesh );//load
            }
            else if ( (x0 >=5   && fabs ( ( yf-500 ) ) <tol)  )
            {
               // new TPZGeoElRefPattern< pzgeom::TPZGeoLinear> ( iel, TopoLine, -5, *gmesh );//toprigth
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
    for ( int d = 0; d<1; d++ )
    {
        int nel = gmesh->NElements();
        TPZManVector<TPZGeoEl *> subels;
        for ( int iel = 0; iel<nel; iel++ )
        {
            TPZGeoEl *gel = gmesh->ElementVec() [iel];
            gel->Divide ( subels );
        }
    }
    std::ofstream print("gmesh.txt");
	gmesh->Print(print);
    std::ofstream files ( "teste-mesh.vtk" );
    TPZVTKGeoMesh::PrintGMeshVTK ( gmesh,files,false );
    cout << "d" << endl;
    return gmesh;

}


TPZCompMesh *Footing::CreateCompMesh(TPZGeoMesh *gmesh) {


    unsigned int dim  = 2;
    const std::string name ( "ElastoPlastic COMP MESH Footing Problem " );

    TPZCompMesh * cmesh =  new TPZCompMesh ( gmesh );
    cmesh->SetName ( name );
    cmesh->SetDefaultOrder (2 );
    cmesh->SetDimModel ( dim );


    //Pu = (2+pi)c = 218.183

      // Mohr Coulomb data
    //REAL mc_cohesion    = 490.;//kPa
    REAL mc_cohesion    = 49.;//N/cm^2
    REAL mc_phi         = 0.0001*M_PI/180;
    REAL mc_psi         = mc_phi;

    /// ElastoPlastic Material using Mohr Coulomb
    // Elastic predictor
    TPZElasticResponse ER;
    REAL nu = 0.48;
    //REAL E = 10000000.;////kPa
    REAL E = 1000000.;////N/cm^2

    LEMC plasticstep;
    ER.SetEngineeringData(E,nu);
    //ER.SetUp ( E, nu );
    plasticstep.fER =ER;
    // LEMC.SetElasticResponse( ER );
    plasticstep.fYC.SetUp ( mc_phi, mc_psi, mc_cohesion, ER );

    plasticstep.fER = ER;
    int PlaneStrain = 1;

    plasticmat * material = new plasticmat ( 1,PlaneStrain );

    REAL factor;
    TPZManVector<REAL, 3> bodyforce ( 3,0. );
    factor=1.;
    bodyforce[1]=0.;

    material->SetPlasticityModel ( plasticstep );

    material->SetId ( 1 );

    //material->SetWhichLoadVector ( 0 ); //option to compute the total internal force vecor fi=(Bt sigma+ N (b+gradu))

   // material->SetLoadFactor ( factor );

    material->SetBodyForce ( bodyforce );

    cmesh->InsertMaterialObject(material);

       TPZFMatrix<STATE>  val1 ( 2,2,0. );

    TPZManVector<STATE,2> val2 ( 2,0. );

    val2[1]=1.;
    auto * bcleft = material->CreateBC(material,-1,3,val1,val2);//bottom

    val2[0]=1.;
    val2[1]=0.;
    auto * bcbottom = material->CreateBC(material,-2,3,val1,val2);//rigth

    val2[0]=1.;
    val2[1]=0;
    auto *  bcrigth = material->CreateBC(material,-3,3,val1,val2);//left

    val2 [0] = 0;
    //val2 ( 1,0 ) = -490 ;//kPa
    val2 [1] = -49 ;//N/cm^2
    auto * bcload = material->CreateBC ( material, -4, 1, val1, val2 );//ramp


	cmesh->InsertMaterialObject(bcleft);

    cmesh->InsertMaterialObject(bcbottom);

    cmesh->InsertMaterialObject(bcrigth);

    cmesh->InsertMaterialObject(bcload);

	cmesh->SetAllCreateFunctionsContinuousWithMem();

    cmesh->AutoBuild();

    std::ofstream print("cmeshf.txt");
	cmesh->Print(print);

    return cmesh;
}

void Footing::PostProcessVariables ( TPZStack<std::string> &scalNames, TPZStack<std::string> &vecNames )
{




    scalNames.Push ( "PlStrainSqJ2" );
    scalNames.Push ( "PlStrainSqJ2El" );
    scalNames.Push ( "PlAlpha" );

    scalNames.Push ( "DisplacementX" );
    scalNames.Push ( "DisplacementY" );
    scalNames.Push ( "DisplacementZ" );
    vecNames.Push ( "Displacement" );
    vecNames.Push ( "PrincipalStress" );


}

void Footing::PostPlasticity(string vtkd,TPZCompMesh * cmesh)
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
void  Footing::CreatePostProcessingMesh (TPZPostProcAnalysis * PostProcess ,TPZCompMesh * cmesh)
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

void Footing::findid(TPZCompMesh *cmesh,TPZVec<REAL> coord) {

    std::ofstream print("loadvsdisplacement.dat");
    TPZGeoMesh * gmesh =  cmesh->Reference();

    int dim =2;

    int id;
    int nels = gmesh->NElements();
    for(int iel=0; iel<nels; iel++) {


        TPZVec<REAL> qsi(dim,0.);
        TPZGeoEl * gel = gmesh->ElementVec()[iel];

        TPZCompEl * cel = gel->Reference();



        if(gel->MaterialId()<0)continue;
        bool check = gel->ComputeXInverse(coord,qsi,1.e-5);

        if(check==true)
        {

            int nnodes=gel->NNodes();
            for(int inode=0; inode<nnodes; inode++)
            {
                TPZVec<REAL> co(dim);
                TPZGeoNode* node = gel->NodePtr(inode);
                node->GetCoordinates(co);
                id = node->Id();

                long elementid1;
                if(fabs(co[0]-coord[0])<1.e-3 && fabs(co[1]-coord[1])<1.e-3)
                {

                }
            }


        }
    }

}

