// slopeconfigure.h
#ifndef SLOPE_H // include guard
#define SLOPE_H

#include <cmath>
#include <set>
#include <iostream>
#include <fstream>
#include <string>
using std::chrono::high_resolution_clock;
using std::chrono::duration_cast;
using std::chrono::duration;
using std::chrono::milliseconds;
using std::chrono::seconds;
#include "pzblockdiag.h"
#include "TPZSpStructMatrix.h"
#include "pzfstrmatrix.h"
#include "pzbdstrmatrix.h"
#include "TPZFileStream.h"
#include <TPZBFileStream.h>
typedef TPZPlasticStepPV<TPZYCMohrCoulombPV, TPZElasticResponse> plasticmorh;
typedef   TPZMatElastoPlastic2D < TPZPlasticStepPV<TPZYCMohrCoulombPV, TPZElasticResponse>, TPZElastoPlasticMem > plasticmat;
//typedef TPZNonLinearAnalysis typedefanal;
typedef TPZElastoPlasticAnalysis typedefanal;
using namespace std;

class Slope
{
public:

    enum  SolverType { EStep, EPardiso};

    //contrutor que inicializa os dados da analise do talude
    Slope( TPZGeoMesh * gmesh,int porder,int ref, REAL gammaagua, REAL gammasolo,REAL coes,REAL atrito);

    //contrutor que inicializa os dados da analise do talude
    TPZCompMesh * CreateCMeshElastoplastic ( TPZGeoMesh *gmesh, int pOrder );
    TPZCompMesh *CreateCompMesh(TPZGeoMesh *gmesh,int porder) ;
    TPZGeoMesh *  CreateGMeshSlope ( int ref );
    typedefanal *  SlopeAnalysis(int type,int numthreads);

    void ApplyGravityLoad(TPZManVector<REAL, 3> bodyforce);

    void PostPlasticity(string vtkd);

    void  CreatePostProcessingMesh (TPZPostProcAnalysis * PostProcess );

    void PostProcessVariables ( TPZStack<std::string> &scalNames, TPZStack<std::string> &vecNames );

    void InitializeSimulation(int porder,int ref);

    REAL Solve(int imc);

    void LoadingRamp (  REAL factor );

    REAL ShearRed ( int maxcout,REAL FS0,REAL fstol );

   // REAL ShearRed (int imc );

   void DivideElementsAbove(REAL refineaboveval, std::set<long> &elindices);

   void PRefineElementsAbove(REAL refineaboveval, int porder, std::set<long> &elindices);

    // Get the vector of element plastic deformations
    void ComputeElementDeformation();

    void ApplyHistory(std::set<long> &elindices);

    int FindElement ( TPZCompMesh * cmesh,TPZManVector<REAL,3>&vecx, TPZManVector<REAL,3>&vecxi );

    void TransferFieldsSolutionFrom ( int isol );

    void ShearReductionIntegrationPoints(REAL FS);

    void InitializeMemory();
    void InitializePOrderInRegion();

    void Write(TPZStream &out);
/// Read the data from the input stream
    void Read(TPZStream &input);

    void SetCompMeshField(TPZVec<TPZCompMesh*>  CompMeshVec)
    {
        fCompMeshVecField =CompMeshVec;
    }

public:
    TPZCompMesh *fCompMesh;
    TPZCompMesh *fCompMeshField;

protected:

    TPZVec<TPZCompMesh*> fCompMeshVecField;
    int fPorder;
    REAL fgammaagua;
    REAL fgammasolo;
    TPZGeoMesh *fGmesh;
    int fref;
    REAL fCoes;
    REAL fAtrito;
    int fSolver;
   // TPZGeoMesh *fGmesh2;
   // typedefanal *fAnalysis;

    TPZVec<REAL> fPlasticDeformSqJ2;


};



Slope::Slope( TPZGeoMesh * gmesh,int porder,int ref,REAL gammaagua, REAL gammasolo,REAL coes, REAL atrito)
{
    fCoes=coes;
    fAtrito=atrito;
    fgammaagua=gammaagua;
    fgammasolo=gammasolo;
    fGmesh = gmesh;
    fPorder=porder;
    fref = ref;
    fCompMesh = CreateCMeshElastoplastic ( fGmesh, fPorder );
    InitializeMemory();
    //InitializePOrderInRegion();
    fSolver=0;

}

REAL Slope::Solve(int imc)
{
    int neqold;
    REAL FSOLD=1000.;
    fCompMesh = CreateCMeshElastoplastic ( fGmesh, fPorder );

    int neq=fCompMesh->NEquations();
    cout << "NUMBER OF EQUATIONS  = " << neq << endl;
    TransferFieldsSolutionFrom(imc);
    InitializePOrderInRegion();
    REAL FS = ShearRed(20,0.5,0.01);

    cout << "Refining.."<<endl;
    std::set<long> elindices,elindices2;
    for (int iref=1; iref<=fref; iref++ )
    {
        ComputeElementDeformation();
        PRefineElementsAbove(0.01, fPorder+iref,elindices2);
        DivideElementsAbove ( 0.01,elindices );
        neqold=neq;
        neq=fCompMesh->NEquations();
        TransferFieldsSolutionFrom(imc);
        cout << "# of equations  = " <<neq << " fabs(FS-FSOLD)  "  << fabs(FS-FSOLD)  << endl;
        FSOLD=FS;
        FS = ShearRed(20,0.5,0.01);

        //if(neq>5000||neqold==neq||fabs(FS-FSOLD)<0.001)
       // {
          //  cout << " # of equations execeded the maximum, or no element was refined, exiting refinement method."<<endl;
         //   break;
        //}
    }

    auto var=to_string ( imc );
    string meshref = "/home/diogo/Dropbox/adaptive-random-fields-applyed-to-slopes/results/result2/post/refinidemesh-grid";
    meshref+=var;
    meshref+=".vtk";
    std::ofstream files ( meshref );
    TPZVTKGeoMesh::PrintGMeshVTK (fCompMesh->Reference(),files,true );

    return FS;
}

void Slope::InitializePOrderInRegion()
{
//     fGmesh->ResetReference();
//     fCompMesh->LoadReferences();
//     long nelem = fCompMesh->NElements();
//     for (long el=0; el<nelem; el++) {
//         TPZCompEl *cel = fCompMesh->ElementVec()[el];
//         TPZInterpolationSpace *intel = dynamic_cast<TPZInterpolationSpace *>(cel);
//         TPZManVector<REAL,3> point ( 2,0. );
//         TPZMaterialDataT<REAL> data;
//         intel->InitMaterialData ( data );
//         intel->ComputeRequiredData ( data, point );
//
//         if(data.x[0]<10)
//         {
//             intel->SetPreferredOrder(5);
//         }else{
//             //intel->SetPreferredOrder(2);
//         }
//     }
//
//     fCompMesh->AdjustBoundaryElements();
//     fCompMesh->InitializeBlock();
}

void Slope::InitializeMemory()
{
     TPZMatWithMem<TPZElastoPlasticMem> *pMatWithMem2 = dynamic_cast<TPZMatWithMem<TPZElastoPlasticMem> *> ( fCompMesh->MaterialVec() [1] );
    if ( pMatWithMem2 )
    {
        pMatWithMem2->SetUpdateMem ( true );
    }
    else
    {
        DebugStop();
    }


    int nels =  fCompMesh->NElements();

    for ( int iel=0; iel<nels; iel++ )
    {

        TPZCompEl *cel = fCompMesh->ElementVec() [iel];
        TPZInterpolationSpace *intel = dynamic_cast<TPZInterpolationSpace *> ( cel );
        if ( !cel || !intel || dynamic_cast<TPZMatWithMem<TPZElastoPlasticMem> *> ( intel->Material() )!= pMatWithMem2 ||(intel->Material()->Id()<1))
        {
            continue;
        }


        const TPZIntPoints &intpoints = intel->GetIntegrationRule();
        int nint = intpoints.NPoints();
        TPZManVector<REAL,3> point ( 2,0. );


        TPZMaterialDataT<REAL> data;
        intel->InitMaterialData ( data );
        data.fNeedsSol = false;

        for ( long ip =0; ip<nint; ip++ )
        {
            REAL weight;
            intpoints.Point ( ip, point, weight );
            data.intLocPtIndex = ip;
            intel->ComputeRequiredData ( data, point );

            int indexplastic =data.intGlobPtIndex;
            TPZElastoPlasticMem &mem = pMatWithMem2->MemItem(indexplastic);
            mem.m_elastoplastic_state.fmatpropinit.Resize ( 3 );
            mem.m_elastoplastic_state.fmatprop.Resize ( 3 );

            mem.m_elastoplastic_state.fmatpropinit[0]=fCoes;
            mem.m_elastoplastic_state.fmatpropinit[1]=fAtrito;
            mem.m_elastoplastic_state.fmatprop[0]=fCoes;
            mem.m_elastoplastic_state.fmatprop[1]=fAtrito;

        }


    }
    pMatWithMem2->SetUpdateMem ( false );
}
void Slope::ShearReductionIntegrationPoints(REAL FS)
{
    TPZMatWithMem<TPZElastoPlasticMem> *pMatWithMem2 = dynamic_cast<TPZMatWithMem<TPZElastoPlasticMem> *> ( fCompMesh->MaterialVec() [1] );
    if ( pMatWithMem2 )
    {
        pMatWithMem2->SetUpdateMem ( true );
    }
    else
    {
        DebugStop();
    }


    int nels =  fCompMesh->NElements();

    for ( int iel=0; iel<nels; iel++ )
    {

        TPZCompEl *cel = fCompMesh->ElementVec() [iel];
        TPZInterpolationSpace *intel = dynamic_cast<TPZInterpolationSpace *> ( cel );
        if ( !cel || !intel || dynamic_cast<TPZMatWithMem<TPZElastoPlasticMem> *> ( intel->Material() )!= pMatWithMem2 ||(intel->Material()->Id()<1))
        {
            continue;
        }


        const TPZIntPoints &intpoints = intel->GetIntegrationRule();
        int nint = intpoints.NPoints();
        TPZManVector<REAL,3> point ( 2,0. );


        TPZMaterialDataT<REAL> data;
        intel->InitMaterialData ( data );
        data.fNeedsSol = true;

        for ( long ip =0; ip<nint; ip++ )
        {
            REAL weight;
            intpoints.Point ( ip, point, weight );
            data.intLocPtIndex = ip;
            intel->ComputeRequiredData ( data, point );

            int indexplastic =data.intGlobPtIndex;
            TPZElastoPlasticMem &mem = pMatWithMem2->MemItem(indexplastic);
            //mem.m_elastoplastic_state.fmatprop.Resize ( 3 );
            if(!mem.m_elastoplastic_state.fmatprop.size())
            {
                cout << "deve-se inicializar corretamente o matprop"<<endl;
                DebugStop();
            }
            REAL coes0=mem.m_elastoplastic_state.fmatpropinit[0];
            REAL atrito0=mem.m_elastoplastic_state.fmatpropinit[1];
            mem.m_elastoplastic_state.fmatprop[0]=coes0/FS;
            mem.m_elastoplastic_state.fmatprop[1]=atan ( tan ( atrito0 ) /FS );

        }


    }
    pMatWithMem2->SetUpdateMem ( false );
}


//Tranfere a solucao nodal da malha cmesh para os pontos de integracao da malha fCompMesh. Este metodo e usado para transferir a solucao dos
//campos estocasticos. O metodo findelement e caro, e custa muito ao monte carlo.
void Slope::TransferFieldsSolutionFrom ( int isol )
{

    TPZMatWithMem<TPZElastoPlasticMem> *pMatWithMem2 = dynamic_cast<TPZMatWithMem<TPZElastoPlasticMem> *> ( fCompMesh->MaterialVec() [1] );
    if ( pMatWithMem2 )
    {
        pMatWithMem2->SetUpdateMem ( true );
    }
    else
    {
        DebugStop();
    }


    int nels =  fCompMesh->NElements();

    int ncmeshes= fCompMeshVecField.size();
    TPZVec<TPZMaterial *> matvec(ncmeshes);
    TPZVec<TPZGeoMesh *> geomeshvec(ncmeshes);
    for(int imesh=0;imesh<ncmeshes;imesh++)
    {
        matvec[imesh] = fCompMeshVecField[imesh]->FindMaterial ( 1 );
        geomeshvec[imesh] =fCompMeshVecField[imesh]->Reference();
    }

    for ( int iel=0; iel<nels; iel++ )
    {

        TPZCompEl *cel = fCompMesh->ElementVec() [iel];
        TPZInterpolationSpace *intel = dynamic_cast<TPZInterpolationSpace *> ( cel );
        if ( !cel || !intel || dynamic_cast<TPZMatWithMem<TPZElastoPlasticMem> *> ( intel->Material() )!= pMatWithMem2 ||(intel->Material()->Id()<1))
        {
            continue;
        }


        const TPZIntPoints &intpoints = intel->GetIntegrationRule();
        int nint = intpoints.NPoints();
        TPZManVector<REAL,3> point ( 2,0. );
       // TPZVec<REAL> point ( 3,0. );


        TPZMaterialDataT<REAL> data;
        intel->InitMaterialData ( data );
        data.fNeedsSol = true;



        for ( long ip =0; ip<nint; ip++ )
        {
            REAL weight;
            intpoints.Point ( ip, point, weight );
            data.intLocPtIndex = ip;
            intel->ComputeRequiredData ( data, point );

            int indexplastic =data.intGlobPtIndex;
            TPZElastoPlasticMem &mem = pMatWithMem2->MemItem(indexplastic);
            mem.m_elastoplastic_state.fmatprop.Resize ( 3 );
            mem.m_elastoplastic_state.fmatpropinit.Resize ( 3 );

            for(int imesh=0;imesh<ncmeshes;imesh++)
            {

                long elementid1 = 0;
                TPZManVector<REAL,3> qsi ( 2,0. );
                //TPZVec<REAL> qsi(3,0.);
                int targetdim=2;
                TPZCompEl*celfield;
                TPZGeoEl *gelfield;
                if(true)
                {
                   elementid1= FindElement ( fCompMeshVecField[imesh],data.x,qsi );
                    celfield =  fCompMeshVecField[imesh]->ElementVec()[elementid1];

                }else{
                    gelfield = geomeshvec[imesh]->FindElement ( data.x, qsi, elementid1,targetdim );
                    celfield = gelfield->Reference();

                }
                TPZInterpolationSpace *intelfield = dynamic_cast<TPZInterpolationSpace *> ( celfield );


                TPZMaterialDataT<REAL> datafield;
                datafield.fNeedsSol = true;
                intelfield->InitMaterialData ( datafield );
                datafield.fNeedsSol = true;
                intelfield->ComputeRequiredData ( datafield, qsi );

                REAL datarealvalue=datafield.sol[isol][0];

                mem.m_elastoplastic_state.fmatpropinit[imesh] = datarealvalue;
                mem.m_elastoplastic_state.fmatprop[imesh] = datarealvalue;
            }
        }
        pMatWithMem2->SetUpdateMem ( false );

    }

}


int Slope::FindElement ( TPZCompMesh * cmesh,TPZManVector<REAL,3>&vecx, TPZManVector<REAL,3>&vecxi )
{
        int nels = cmesh->NElements();

        REAL tol=1.e-5;

        for ( int iel=0; iel<nels; iel++ ) {

                TPZCompEl * cel =  cmesh->ElementVec() [iel];
                if ( !cel ) {
                        continue;
                }
                TPZGeoEl * gel = cel->Reference();
                if ( !gel ) {
                        continue;
                }
                bool find = gel->ComputeXInverse ( vecx,vecxi,tol );
                if ( find ) {
                        //cout << "ponto encontrado no elemento id "<<iel <<endl;
                        //cout << "coordenadas no elemento mestre = "<< vecxi << endl;
                        return iel;
                } else {
                        //cout << " ponto NAO encontrado! " <<endl;
                        //DebugStop();
                }

        }

            cout << " int KLAnalysis::FindElement  ponto NAO encontrado! " <<endl;
            DebugStop();

}


TPZGeoMesh *  Slope::CreateGMeshSlope ( int ref )
{
    const std::string name ( "Darcy Flow Slope" );

    TPZGeoMesh *gmesh  =  new TPZGeoMesh();

    gmesh->SetName ( name );

    gmesh->SetDimension ( 2 );

    TPZVec<REAL> coord ( 2 );

vector<vector<double>> co= {
    {0., 0.}, {75., 0.}, {75., 30.},{45., 30.},
    {35., 40.}, {0.,40.},{35./3., 40.},{2 * 35/3., 40.},
    {30., 40.},{30., 30.}, {60.,30.},{2* 35./3.,2* 35/3.},
    {45., 2* 35/3.},{35./3., 35/3.}, {60., 35./3.}
};
vector<vector<int>> topol = {
    {0,  1,  14, 13},{1,  2,  10, 14}, {14, 10, 3,  12},
    {13, 14, 12, 11},{11, 12, 3,  9}, {9,  3,  4,  8},
    {11, 9,  8,  7},{13, 11, 7, 6},{0, 13,  6, 5}
};

    gmesh->NodeVec().Resize ( co.size() );

    for ( int inode=0; inode<co.size(); inode++ ) {
        coord[0] = co[inode][0];
        coord[1] = co[inode][1];
        gmesh->NodeVec() [inode] = TPZGeoNode ( inode, coord, *gmesh );
    }
    TPZVec <long> TopoQuad ( 4 );
    for ( int iel=0; iel<topol.size(); iel++ ) {
        TopoQuad[0] = topol[iel][0];
        TopoQuad[1] = topol[iel][1];
        TopoQuad[2] =	topol[iel][2];
        TopoQuad[3] = topol[iel][3];
        new TPZGeoElRefPattern< pzgeom::TPZGeoQuad> ( iel, TopoQuad, 1,*gmesh );
    }




    int id = topol.size();
    TPZVec <long> TopoLine ( 2 );
    TopoLine[0] = 0;
    TopoLine[1] = 1;
    new TPZGeoElRefPattern< pzgeom::TPZGeoLinear> ( id, TopoLine, -1, *gmesh );//bottom

    id++;
    TopoLine[0] = 1;
    TopoLine[1] = 2;
    new TPZGeoElRefPattern< pzgeom::TPZGeoLinear> ( id, TopoLine, -2, *gmesh );//rigth

    id++;
    TopoLine[0] = 2;
    TopoLine[1] = 3;
    new TPZGeoElRefPattern< pzgeom::TPZGeoLinear> ( id, TopoLine, -3, *gmesh );//top-rigth

    id++;
    TopoLine[0] = 4;
    TopoLine[1] = 5;
    new TPZGeoElRefPattern< pzgeom::TPZGeoLinear> ( id, TopoLine, -4, *gmesh ); //top-left

    id++;
    TopoLine[0] = 3;
    TopoLine[1] = 4;
    new TPZGeoElRefPattern< pzgeom::TPZGeoLinear> ( id, TopoLine, -6, *gmesh ); //ramp

    id++;
    TopoLine[0] = 5;
    TopoLine[1] = 0;
    new TPZGeoElRefPattern< pzgeom::TPZGeoLinear> ( id, TopoLine, -5, *gmesh ); //left





    gmesh->BuildConnectivity();
    for ( int d = 0; d<ref; d++ ) {
        int nel = gmesh->NElements();
        TPZManVector<TPZGeoEl *> subels;
        for ( int iel = 0; iel<nel; iel++ ) {
            TPZGeoEl *gel = gmesh->ElementVec() [iel];
            gel->Divide ( subels );
        }
    }

    string meshref = "gmesh.vtk";
    std::ofstream files ( meshref );
    TPZVTKGeoMesh::PrintGMeshVTK ( gmesh,files,true );
    return gmesh;
}

void Slope::LoadingRamp (  REAL factor )
{
    plasticmat * body= dynamic_cast<plasticmat *> ( fCompMesh->FindMaterial ( 1 ) );
    //plasticmatcrisfield * body= dynamic_cast<plasticmatcrisfield *> ( cmesh->FindMaterial ( 1 ) );
    TPZManVector<REAL, 3> force ( 3,0. );

    force[1]= ( fgammaagua-fgammasolo );

    //force[1]=(-gammasolo);
    body->SetLoadFactor ( factor );
    body->SetBodyForce ( force );

}


REAL Slope::ShearRed (int maxcout,REAL FS0,REAL fstol )
{
    LoadingRamp(1.);

    REAL FS=FS0,FSmax=10.,FSmin=0.,tol=fstol;
    int neq = fCompMesh->NEquations();
    int counterout = 0;
    int numthreads = 10;
    int type =fSolver;
    bool conv=false;
    auto t0 =high_resolution_clock::now();
    REAL FSN=1000;

    do {

        fCompMesh->Solution().Zero();
        typedefanal* anal =  SlopeAnalysis(type,numthreads);
        REAL norm = 1000.;
        REAL tol2 = 1.e-3;
        int NumIter = 100;
        bool linesearch = true;
        bool checkconv = false;
        int iters;


        //std::cout << "start solving iterative process"<< endl ;

        ShearReductionIntegrationPoints(FS);
        auto t1 = high_resolution_clock::now();
        conv  =anal->IterativeProcess2(cout,tol2, NumIter,  linesearch,  checkconv,iters);
        //conv  =anal->FindRoot(iters);

        auto t2 = high_resolution_clock::now();
        auto ms_int = duration_cast<milliseconds> ( t2 - t1 );
        norm = Norm ( anal->Rhs() );
        cout << "| step = " << counterout << " FS = "<< FS <<" tempo  iterproc = "<<ms_int.count() << " ms " << " conv?" << conv << " iters = " <<iters<< endl;
         //conv = anal->FindRoot ( );


        FSN=FS;
        if ( conv==false ) {

            FSmax = FS;
            FS = ( FSmin + FSmax ) / 2.;
        } else {

            FSmin = FS;
			FS = 1. / ( ( 1. / FSmin + 1. / FSmax ) / 2. );
        }
        if(fabs(FSN-FS)<1.e-3 && conv==true)
        {
            anal->AcceptSolution();
            //conv=true;
        }

        counterout++;
        if(( FSmax - FSmin ) / FS < tol  && conv==true)
        {
            anal->AcceptSolution();
            //conv=true;
        }
        delete anal;
    }  while ( (( FSmax - FSmin ) / FS > tol || conv==false) && counterout<maxcout);

        auto t3 = high_resolution_clock::now();
        auto timeinmili = duration_cast<seconds> ( t3 - t0 );
        std::cout << "final safety factor "<< FS << " total time in ShearRed = "<< timeinmili.count() << " s "<<std::endl;
        return ( FSmax + FSmin )/2;
}

/*
REAL Slope::ShearRed (int imc )
{
    int numthreads = 10;
    int type =EStep;
    //std::cout << "aqqqqq "<<std::endl;
    LoadingRamp(1.);

    REAL FS=0.8,FSmax=5.,FSmin=0.,tol=0.001;
    int neq = fCompMesh->NEquations();

    TPZFMatrix<REAL> displace(neq,1),displace0(neq,1);

    int counterout = 0;

    plasticmat *material = dynamic_cast<plasticmat *> ( fCompMesh->MaterialVec() [1] );
    TPZPlasticStepPV<TPZYCMohrCoulombPV, TPZElasticResponse> LEMC = material->GetPlasticity();
    TPZElasticResponse ER = LEMC.fER;
    REAL FSN=1000;
    bool conv=false;
    auto t0 =high_resolution_clock::now();
    do {

        fCompMesh->Solution().Zero();
        typedefanal* anal =  SlopeAnalysis(type,numthreads);

        //TransferSolutionFromShearRed (  fCompMeshField,imc, FS );
        auto t1 = high_resolution_clock::now();

        //TransferSolutionFrom (fCompMeshField,imc,FS );

        auto t2 = high_resolution_clock::now();
        auto ms_int = duration_cast<milliseconds> ( t2 - t1 );

        cout <<" tempo  transfering = "<<ms_int.count() << endl;

        REAL norm = 1000.;
        REAL tol2 = 0.01;
        int NumIter = 50;
        bool linesearch = true;
        bool checkconv = false;
        int iters;



        t1 = high_resolution_clock::now();
        conv  =anal->IterativeProcess2(cout,tol2, NumIter,  linesearch,  checkconv,iters);
        //conv  =anal->FindRoot(iters);

        t2 = high_resolution_clock::now();
        ms_int = duration_cast<milliseconds> ( t2 - t1 );
        //std::cout << "tempo total iterative process = "<<ms_int.count() << " s\n";
        norm = Norm ( anal->Rhs() );
        //conv  = anal->IterativeProcess ( cout, tol2, NumIter,linesearch,checkconv,iters);
cout << "| step = " << counterout << " FS = "<< FS <<" tempo  iterproc = "<<ms_int.count() << " ms " << " conv?" << conv << " iters = " <<iters<< endl;



        FSN=FS;
        if ( conv==false ) {

            FSmax = FS;
            FS = ( FSmin + FSmax ) / 2.;
        } else {

            FSmin = FS;
			FS = 1. / ( ( 1. / FSmin + 1. / FSmax ) / 2. );
        }
        if(fabs(FSN-FS)<1.e-3)
        {
            anal->AcceptSolution();
            conv=true;
        }
        counterout++;
        if(( FSmax - FSmin ) / FS < tol)
        {
            anal->AcceptSolution();
        }
    }  while ( ( FSmax - FSmin ) / FS > tol || conv==false);

        auto t3 = high_resolution_clock::now();
        auto timeinmili = duration_cast<seconds> ( t3 - t0 );
        std::cout << "final safety factor "<< FS << " total time in ShearRed = "<< timeinmili.count() << " s "<<std::endl;
        return ( FSmax + FSmin )/2;
}
*/



TPZCompMesh *Slope::CreateCompMesh(TPZGeoMesh *gmesh,int porder) {


    unsigned int dim  = porder;
    const std::string name ( "ElastoPlastic COMP MESH Slope Problem " );

    TPZCompMesh * cmesh =  new TPZCompMesh ( gmesh );
    cmesh->SetName ( name );
    cmesh->SetDefaultOrder (porder );
    cmesh->SetDimModel ( dim );


    //Pu = (2+pi)c = 218.183

      // Mohr Coulomb data
    //REAL mc_cohesion    = 490.;//kPa
    REAL mc_cohesion    = fCoes;//N/cm^2
    REAL mc_phi         = fAtrito;
    REAL mc_psi         = mc_phi;

    /// ElastoPlastic Material using Mohr Coulomb
    // Elastic predictor
    TPZElasticResponse ER;
    REAL nu = 0.49;
    //REAL E = 10000000.;////kPa
    REAL E = 20000.;////N/cm^2

    plasticmorh plasticstep;
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

    material->SetWhichLoadVector ( 0 ); //option to compute the total internal force vecor fi=(Bt sigma+ N (b+gradu))

    material->SetLoadFactor ( factor );

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
    auto *  bcrigth = material->CreateBC(material,-5,3,val1,val2);//left


	cmesh->InsertMaterialObject(bcleft);

    cmesh->InsertMaterialObject(bcbottom);

    cmesh->InsertMaterialObject(bcrigth);



	cmesh->SetAllCreateFunctionsContinuousWithMem();

    cmesh->AutoBuild();

    std::ofstream print("cmeshslopecphi.txt");
	cmesh->Print(print);

    return cmesh;
}

TPZCompMesh * Slope::CreateCMeshElastoplastic ( TPZGeoMesh *gmesh, int pOrder )
{
	    // Creating computational mesh:
    TPZCompMesh * cmesh = new TPZCompMesh ( gmesh );

    //cmesh->Print(std::cout);
	cmesh->SetDefaultOrder ( pOrder );
    //cmesh->GetDefaultOrder

	int dim = 2 ;

	cmesh->SetDimModel ( dim );



	int matid=1;
	int planestrain=1;

	//auto * material = new plasticmat(matid,planestrain);

    REAL E=20000;
    REAL nu=0.49;
    TPZElasticResponse  elasticresponse;

    //REAL lambda = nu * E / ((1. + nu)*(1. - 2. * nu));
    //REAL mu = E / (2. * (1. + nu));
    elasticresponse.SetEngineeringData(E,nu);

    // Mohr Coulomb data
    REAL mc_cohesion    = fCoes;//kPa
    REAL mc_phi         = fAtrito;
    REAL mc_psi         = mc_phi;

    //elasticresponse.Print(std::cout);

    plasticmorh mohrcoulombplasticstep;

    mohrcoulombplasticstep.fYC.SetUp ( mc_phi, mc_psi, mc_cohesion, elasticresponse );

    mohrcoulombplasticstep.fER = elasticresponse;

    //mohrcoulombplasticstep.Print(std::cout);

    int PlaneStrain = 1;

    plasticmat * material = new plasticmat ( matid,PlaneStrain );

    material->SetPlasticityModel(mohrcoulombplasticstep);

    material->SetId(matid);

    REAL factor;
    TPZManVector<REAL, 3> bodyforce ( 3,0. );
    factor=1.;
    bodyforce[1]=-20.;

   // material->SetPlasticity ( plasticstep );

    material->SetId ( 1 );

    material->SetWhichLoadVector ( 0 ); //option to compute the total internal force vecor fi=(Bt sigma+ N (b+gradu))

    material->SetLoadFactor ( factor );

    material->SetBodyForce ( bodyforce );
   // material->

    //material->Print(std::cout);

    cmesh->InsertMaterialObject ( material );

    //cmesh->Print(std::cout);

    // boundary condition
    TPZFMatrix<STATE>  val1 ( 2,2,0. );

    TPZManVector<STATE,2> val2 ( 2,0. );

	int directionaldirichlet = 3 ;

	val2[0]=1;
	val2[1]=1;
    auto * BCond0 = material->CreateBC ( material, -1, directionaldirichlet, val1, val2 );

	val2[0]=1;
	val2[1]=0;
    auto * BCond1 = material->CreateBC ( material, -2, directionaldirichlet, val1, val2 );

	val2[0]=1;
	val2[1]=0;
	auto * BCond2 = material->CreateBC ( material, -5, directionaldirichlet, val1, val2 );

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


void Slope::ApplyGravityLoad(TPZManVector<REAL, 3> bodyforce)
{
    plasticmat * body= dynamic_cast<plasticmat *> ( fCompMesh->FindMaterial ( 1 ) );
    body->SetBodyForce ( bodyforce );

}

typedefanal *  Slope::SlopeAnalysis(int type,int numthreads)
{

    typedefanal *analysis = new typedefanal(fCompMesh,cout);
   // analysis->SetCompMesh(fCompMesh,true);
switch(type) {
    case 0:
        {
			//cout << "Solver called with TPZStepSolver\n";
            TPZSkylineStructMatrix<STATE> matskl(fCompMesh);
            matskl.SetNumThreads ( numthreads );
             analysis->SetStructuralMatrix ( matskl );
             TPZStepSolver<STATE> step;
             step.SetDirect ( ELDLt );
             analysis->SetSolver ( step );
             break;
        }
    case 1:
        {
            //cout << "Solver called with TPZPardisoSolver\n";
            TPZSSpStructMatrix<STATE> SSpStructMatrix(fCompMesh);
            SSpStructMatrix.SetNumThreads ( numthreads );
            analysis->SetStructuralMatrix(SSpStructMatrix);
            TPZPardisoSolver<REAL> *pardiso = new TPZPardisoSolver<REAL>;
            analysis->SetSolver(*pardiso);
            break;
        }
    case 2:
    {
        int numiter=1000;
        REAL tol =1.e-6;
        TPZSpStructMatrix<STATE> StrMatrix(fCompMesh);
        analysis->SetStructuralMatrix(StrMatrix);
        StrMatrix.SetNumThreads ( numthreads );
        TPZMatrix<REAL> * mat = StrMatrix.Create();

        TPZBlockDiagonalStructMatrix<STATE> strBlockDiag(fCompMesh);
        TPZStepSolver<REAL> Pre;
        TPZBlockDiagonal<REAL> * block = new TPZBlockDiagonal<REAL>();

        strBlockDiag.AssembleBlockDiagonal(*block); // just to initialize structure
        Pre.SetMatrix(block);
        Pre.SetDirect(ELU);
        TPZStepSolver<REAL> Solver;
        Solver.SetBiCGStab(numiter, Pre, tol, 0);

        Solver.SetMatrix(mat);
        analysis->SetSolver(Solver);
        analysis->SetPrecond(Pre);
        break;
    }
        case 3:
    {
                int numiter=1000;
        REAL tol =1.e-3;
	TPZSpStructMatrix<STATE> StrMatrix(fCompMesh);
    StrMatrix.SetNumThreads ( numthreads );
    analysis->SetStructuralMatrix(StrMatrix);
	TPZMatrix<REAL> * mat = StrMatrix.Create();

    TPZBlockDiagonalStructMatrix<STATE> strBlockDiag(fCompMesh);
    TPZStepSolver<REAL> Pre;
    TPZBlockDiagonal<REAL> * block = new TPZBlockDiagonal<REAL>();

    strBlockDiag.AssembleBlockDiagonal(*block); // just to initialize structure
	Pre.SetMatrix(block);
    Pre.SetDirect(ELU);
    TPZStepSolver<REAL> Solver;
 	Solver.SetBiCGStab(numiter, Pre, tol, 0);
    Solver.SetMatrix(mat);
    analysis->SetSolver(Solver);
	analysis->SetPrecond(Pre);
        break;
    }
        case 4:
        {
			//cout << "Solver called with TPZStepSolver\n";
            TPZSpStructMatrix<STATE> SSpStructMatrix(fCompMesh);
            SSpStructMatrix.SetNumThreads ( numthreads );
            analysis->SetStructuralMatrix(SSpStructMatrix);
             TPZStepSolver<STATE> step;
             step.SetDirect ( ELU );
             analysis->SetSolver ( step );
             break;
        }
            case 5:
        {
			//cout << "Solver called with TPZStepSolver\n";
            TPZSkylineStructMatrix<STATE> matskl(fCompMesh);
            matskl.SetNumThreads ( numthreads );
             analysis->SetStructuralMatrix ( matskl );
             TPZStepSolver<STATE> step;
             step.SetDirect ( ELDLt );
             analysis->SetSolver ( step );
             break;
        }
    default:
        {
            cout << "Solver was not initialized properly\n";
            DebugStop();
        }
    }
    return analysis;
}



void Slope::PostPlasticity(string vtkd)
{

   TPZPostProcAnalysis * postprocdeter = new TPZPostProcAnalysis();
    CreatePostProcessingMesh ( postprocdeter);


    TPZVec<int> PostProcMatIds ( 1,1 );

    TPZStack<std::string> PostProcVars, scalNames, vecNames;

    PostProcessVariables ( scalNames, vecNames );

    //string vtkd = "postprocessdeter.vtk";
    postprocdeter->DefineGraphMesh ( 2,scalNames,vecNames,vtkd );

    postprocdeter->PostProcess ( 0 );

    delete postprocdeter;
}

void  Slope::CreatePostProcessingMesh (TPZPostProcAnalysis * PostProcess )
{
    if ( PostProcess->ReferenceCompMesh() != fCompMesh )
    {

        PostProcess->SetCompMesh ( fCompMesh );

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
        TPZFStructMatrix<REAL> structmatrix ( PostProcess->Mesh() );
        PostProcess->SetStructuralMatrix ( structmatrix );
        PostProcess->SetPostProcessVariables ( PostProcMatIds, PostProcVars );


    }
    //
    //Chamar com o analysis e nao com o postanalysis pois tem o acumulo de sols
    PostProcess->TransferSolution();

}

void Slope::PostProcessVariables ( TPZStack<std::string> &scalNames, TPZStack<std::string> &vecNames )
{

    scalNames.Push ( "POrder" );
    scalNames.Push ( "Atrito" );
    scalNames.Push ( "Coesion" );
    scalNames.Push ( "StrainPlasticJ2" );
     //scalNames.Push ( "VolHardening" );
     vecNames.Push ( "Displacement" );
     //vecNames.Push ( "ShearPlasticDeformation" );
     //vecNames.Push ( "PlasticDeformation" );


}

/// Divide the element using the plastic deformation as threshold
void Slope::DivideElementsAbove(REAL refineaboveval, std::set<long> &elindices)
{

    //int porder =fPorder+3;
    //fGmesh->ResetReference();
    //fCompMesh->LoadReferences();
    TPZManVector<REAL,3> findel(3,0.),qsi(2,0.);


    long nelem = fCompMesh->NElements();
    for (long el=0; el<nelem; el++) {
        TPZCompEl *cel = fCompMesh->ElementVec()[el];
        if (!cel) {
            continue;
        }

        TPZInterpolationSpace *intel = dynamic_cast<TPZInterpolationSpace *>(cel);
        if (!intel) {
            DebugStop();
        }
        TPZMatWithMem<TPZElastoPlasticMem> *pMatWithMem2 = dynamic_cast<TPZMatWithMem<TPZElastoPlasticMem> *> (cel->Material());
        if (!pMatWithMem2) {
            continue;
        }
        const TPZFMatrix<STATE> &elsol = fCompMesh->ElementSolution();
        if (elsol.Get(el,0) <=refineaboveval) {
            continue;
        }
        //intel->PRefine(3);
        int porder = intel->GetPreferredOrder();
        TPZStack<long> subels;
        long index = cel->Index();


        intel->Divide(index, subels,0);
        for (int is=0; is<subels.size(); is++) {
            elindices.insert(subels[is]);
            TPZCompEl *subcel = fCompMesh->ElementVec()[subels[is]];

            TPZInterpolationSpace *subintel = dynamic_cast<TPZInterpolationSpace *>(subcel);
            if (!subintel) {
                DebugStop();
            }
            subintel->SetPreferredOrder(porder);
        }
    }
    // divide elements with more than one level difference
    bool changed = true;
    while (changed) {
        changed = false;
        std::set<long> eltodivide;
        long nelem = fCompMesh->NElements();
        for (long el=0; el<nelem; el++) {
            TPZCompEl *cel = fCompMesh->ElementVec()[el];
            if (!cel) {
                continue;
            }
            TPZInterpolationSpace *intel = dynamic_cast<TPZInterpolationSpace *>(cel);
            if (!intel) {
                DebugStop();
            }
            TPZGeoEl *gel = cel->Reference();
            if (!gel) {
                DebugStop();
            }
            int ns = gel->NSides();
            for (int is=0; is<ns; is++) {
                TPZGeoElSide gelside(gel, is);
                if (gelside.Dimension() != 1) {
                    continue;
                }
                TPZCompElSide big = gelside.LowerLevelCompElementList2(1);
                if (!big) {
                    continue;
                }
                TPZGeoElSide geobig(big.Reference());
                // boundary elements will be refined by AdjustBoundaryElements
                if (geobig.Element()->Dimension() != 2) {
                    continue;
                }
                if (gel->Level()-geobig.Element()->Level() > 1) {
                    eltodivide.insert(big.Element()->Index());
                }
            }
        }
        std::set<long>::iterator it;
        for (it = eltodivide.begin(); it != eltodivide.end(); it++) {
            changed = true;
            long el = *it;
            TPZCompEl *cel = fCompMesh->ElementVec()[el];
            if (!cel) {
                continue;
            }
            TPZInterpolationSpace *intel = dynamic_cast<TPZInterpolationSpace *>(cel);
            if (!intel) {
                DebugStop();
            }
            TPZMatWithMem<TPZElastoPlasticMem> *pMatWithMem2 = dynamic_cast<TPZMatWithMem<TPZElastoPlasticMem> *> (cel->Material());
            if (!pMatWithMem2) {
                continue;
            }

            int porder = intel->GetPreferredOrder();
            TPZStack<long> subels;
            long index = cel->Index();
            intel->Divide(index, subels,0);
            for (int is=0; is<subels.size(); is++) {
                elindices.insert(subels[is]);
                TPZCompEl *subcel = fCompMesh->ElementVec()[subels[is]];
                TPZInterpolationSpace *subintel = dynamic_cast<TPZInterpolationSpace *>(subcel);
                if (!subintel) {
                    DebugStop();
                }
                subintel->SetPreferredOrder(porder);
            }
        }
    }

//     //ApplyHistory(elindices);
//     ComputeElementDeformation();
//     fCompMesh->AdjustBoundaryElements();
//     fcmesh->InitializeBlock();
//     fCompMesh->Solution().Zero();
//    // fneq=fcmesh->NEquations();
//     fCompMesh->Solution().Resize(0, 0);
//     fCompMesh->Solution().Redim(fCompMesh->NEquations(), 1);
//    // fcmesh->LoadReferences();

    fCompMesh->AdjustBoundaryElements();
    fCompMesh->InitializeBlock();
    TPZMatWithMem<TPZElastoPlasticMem> *pMatWithMem2 = dynamic_cast<TPZMatWithMem<TPZElastoPlasticMem> *> ( fCompMesh->MaterialVec() [1] );
    pMatWithMem2->ResetMemory();
//     fCompMesh->Solution().Zero();
//     fCompMesh->Solution().Resize(0, 0);
//     fCompMesh->Solution().Redim(fCompMesh->NEquations(), 1);


}

/// Change the polynomial order of element using the plastic deformation as threshold
void Slope::PRefineElementsAbove(REAL refineaboveval, int porder, std::set<long> &elindices)
{
    fGmesh->ResetReference();
    fCompMesh->LoadReferences();
    long nelem = fCompMesh->NElements();
    for (long el=0; el<nelem; el++) {
        TPZCompEl *cel = fCompMesh->ElementVec()[el];
        if (!cel) {
            continue;
        }
        TPZInterpolationSpace *intel = dynamic_cast<TPZInterpolationSpace *>(cel);
        if (!intel) {
            DebugStop();
        }
        TPZMatWithMem<TPZElastoPlasticMem> *pMatWithMem2 = dynamic_cast<TPZMatWithMem<TPZElastoPlasticMem> *> (cel->Material());
        if (!pMatWithMem2) {
            continue;
        }

        const TPZFMatrix<STATE> &elsol = fCompMesh->ElementSolution();
        if (elsol(el,0) < refineaboveval) {
            continue;
        }
        //cout << "porder = " << porder << endl;
        TPZStack<long> subels;
        long index = cel->Index();
        elindices.insert(index);
        intel->SetPreferredOrder(porder);
    }

    fCompMesh->AdjustBoundaryElements();
    fCompMesh->InitializeBlock();
//     fCompMesh->Solution().Zero();
//     fCompMesh->Solution().Resize(0, 0);
//     fCompMesh->Solution().Redim(fCompMesh->NEquations(), 1);
}
// Get the vector of element plastic deformations
void Slope::ComputeElementDeformation()
{
    long nelem = fCompMesh->NElements();
    fPlasticDeformSqJ2.resize(nelem);
    fPlasticDeformSqJ2.Fill(0.);
    fCompMesh->ElementSolution().Redim(nelem, 1);
    TPZMatWithMem<TPZElastoPlasticMem> *pMatWithMem2 = dynamic_cast<TPZMatWithMem<TPZElastoPlasticMem> *> (fCompMesh->MaterialVec()[1]);
    plasticmat * elastoplasticmat= dynamic_cast<plasticmat *> ( fCompMesh->FindMaterial ( 1 ) );
    if (!pMatWithMem2) {
        fPlasticDeformSqJ2.Fill(0.);
    }
    else
    {
        for (long el = 0; el<nelem; el++) {
            TPZCompEl *cel = fCompMesh->ElementVec()[el];
            fPlasticDeformSqJ2[el] = 0.;
            if (!cel) {
                continue;
            }
            TPZManVector<long> memindices;
            cel->GetMemoryIndices(memindices);
            int numind = memindices.size();
            REAL sqj2el = 0.;
            REAL phivalplane=0.;
            for (int ind=0; ind<numind; ind++)
            {
                int memoryindex = memindices[ind];
                if (memoryindex < 0) {
                    continue;
                }
                TPZElastoPlasticMem &mem = pMatWithMem2->MemItem(memindices[ind]);
                TPZTensor<REAL> plastic =mem.m_elastoplastic_state.EpsP();
                TPZTensor<REAL> total =mem.m_elastoplastic_state.EpsT();
                TPZVec<REAL> phi;
                elastoplasticmat->GetPlasticity().Phi(total,phi);
                 REAL J2 = plastic.J2();
                 REAL sqj2 = sqrt(J2);
                //REAL val=mem.m_elastoplastic_state.VolHardening();
                sqj2el = max(sqj2,sqj2el);
                phivalplane=phi[0];

            }
            fPlasticDeformSqJ2[el] = sqj2el;
        }
    }
    fCompMesh->SetElementSolution(0, fPlasticDeformSqJ2);
}


void Slope::ApplyHistory(std::set<long> &elindices)
{

//     std::set<long>::iterator it;
//     for (it=elindices.begin(); it != elindices.end(); it++) {
//         long elindex = *it;
//         TPZCompEl *cel = fCurrentConfig.fCMesh.ElementVec()[elindex];
//         TPZMatWithMem<TPZElastoPlasticMem> *pMatWithMem2 = dynamic_cast<TPZMatWithMem<TPZElastoPlasticMem> *> (cel->Material());
//         if (!pMatWithMem2) {
//             DebugStop();
//         }
//         // Reset the memory of the integration points of the element
//         TPZManVector<long> pointindices;
//         cel->GetMemoryIndices(pointindices);
//         long npoints = pointindices.size();
//         for (long ip = 0; ip<npoints; ip++) {
//             long ind = pointindices[ip];
//             pMatWithMem2->ResetMemItem(ind);
//         }
//
//         int confindex = 0;
//         //ApplyDeformation(cel);
//     }
}

/// Write the data to the output stream
void Slope::Write(TPZStream &out)
{
    //fGmesh->Write(out, 0);
    fCompMesh->Write(out, 0);
    out.Write(fgammaagua);
    out.Write(fgammasolo);
    int verify = 83562;
    out.Write(&verify);
}

/// Read the data from the input stream
void Slope::Read(TPZStream &input)
{
    //fGmesh->Read(input, 0);
    //fCompMesh->Read(input, fGmesh);
    input.Read(&fgammaagua);
    input.Read(&fgammasolo);
    int verify = 0;
    input.Read(&verify);
    if (verify != 83562)
    {
        DebugStop();
    }
}
#endif /* SLOPE_H */
