// SPDX-FileCopyrightText: 2024 <copyright holder> <email>
// SPDX-License-Identifier: Apache-2.0

#include "SlopeAnalysis.h"

SlopeAnalysis::SlopeAnalysis()
{

}

SlopeAnalysis::SlopeAnalysis ( const SlopeAnalysis& other )
{
        fCohesion=other.fCohesion;
        fAtrito=other.fAtrito;
        //fPorder=other.fPorder;
        fGammaW=other.fGammaW;
        fGammaS=other.fGammaS;
        fGMesh=other.fGMesh;
        fCompMesh=other.fCompMesh;
}
SlopeAnalysis::SlopeAnalysis ( REAL gammaagua, REAL gammasolo,REAL coes,REAL atrito,int ref0, int porder )
{
        fGMesh = TriGMesh(ref0);
        fCompMesh = CreateCMesh(fGMesh,porder,coes,atrito);
        fCohesion=coes;
        fAtrito=atrito;
        //fPorder=porder;
        fGammaW=gammaagua;
        fGammaS=gammasolo;
}
SlopeAnalysis::~SlopeAnalysis()
{

}
REAL SlopeAnalysis::SolveSingleField ( int ifield )
{
        REAL FSOLD;
        int neqold;
        int neq=fCompMesh->NEquations();
        cout << "NUMBER OF EQUATIONS  = " << neq << endl;
        TransferFieldsSolutionFrom ( ifield );

        REAL FS = ShearRed ( 20,0.5,0.01 );

        cout << "Refining.."<<endl;
        std::set<long> elindices,elindices2;
        for ( int iref=1; iref<=2; iref++ ) {
                ComputeElementDeformation();
                PRefineElementsAbove ( 0.01, fCompMesh->GetDefaultOrder()+iref,elindices2 );
                DivideElementsAbove ( 0.01,elindices );
                neqold=neq;
                neq=fCompMesh->NEquations();
                TransferFieldsSolutionFrom ( ifield );
                cout << "# of equations  = " <<neq << " fabs(FS-FSOLD)  "  << fabs ( FS-FSOLD )  << endl;
                FSOLD=FS;
                FS = ShearRed ( 20,0.5,0.01 );

                if ( fabs ( FS-FSOLD ) <0.01 ) {
                        //cout << " FS-FSOLD = "<< fabs ( FS-FSOLD ) <<endl;
                        //break;
                } else if ( neq>10000||neqold==neq ) {
                        //cout << " neq>10000 = "<< neq <<" neqold = "<< neqold <<endl;
                        //break;
                }
        }

        TransferFieldsSolutionFrom ( ifield );
        auto var=to_string ( ifield );
        string meshref = "post/refinidemesh-grid";
        meshref+=var;
        meshref+=".vtk";
        std::ofstream files ( meshref );
        TPZVTKGeoMesh::PrintGMeshVTK ( fCompMesh->Reference(),files,true );

        return FS;
}

void SlopeAnalysis::ApplyGravityLoad ( TPZManVector<REAL, 3> bodyforce )
{
        plasticmat * body= dynamic_cast<plasticmat *> ( fCompMesh->FindMaterial ( 1 ) );
        body->SetBodyForce ( bodyforce );

}

void SlopeAnalysis::LoadingRamp ( REAL factor )
{
        plasticmat * body= dynamic_cast<plasticmat *> ( fCompMesh->FindMaterial ( 1 ) );
        //plasticmatcrisfield * body= dynamic_cast<plasticmatcrisfield *> ( cmesh->FindMaterial ( 1 ) );
        TPZManVector<REAL, 3> force ( 3,0. );

        force[1]= ( fGammaW-fGammaS );

        //force[1]=(-gammasolo);
        body->SetLoadFactor ( factor );
        body->SetBodyForce ( force );

}


REAL SlopeAnalysis::ShearRed ( int maxcout,REAL FS0,REAL fstol )
{
        LoadingRamp ( 1. );


        REAL FS=FS0,FSmax=10.,FSmin=0.,tol=fstol;
        int counterout = 0;
        bool conv=false;
        auto t0 =chrono::high_resolution_clock::now();
        REAL FSN=1000;
        int type=0;
        int numthreads=15;

        do {

                TPZElastoPlasticAnalysis anal =  SetSlopeAnalysis(type,numthreads);
                fCompMesh->Solution().Zero();
                REAL norm = 1000.;
                REAL tol2 = 1.e-3;
                int NumIter = 100;
                bool linesearch = true;
                bool checkconv = false;
                int iters;

                ShearReductionIntegrationPoints ( FS );

                auto t1 = chrono::high_resolution_clock::now();
                conv  =anal.IterativeProcess2 ( cout,tol2, NumIter,  linesearch,  checkconv,iters );


                auto t2 = chrono::high_resolution_clock::now();
                auto ms_int = chrono::duration_cast<chrono::milliseconds> ( t2 - t1 );
                norm = Norm ( anal.Rhs() );
                cout << "| step = " << counterout << " FS = "<< FS <<" tempo  iterproc = "<<ms_int.count() << " ms " << " conv?" << conv << " iters = " <<iters<< endl;


                FSN=FS;
                if ( conv==false ) {

                        FSmax = FS;
                        FS = ( FSmin + FSmax ) / 2.;
                } else {

                        FSmin = FS;
                        FS = 1. / ( ( 1. / FSmin + 1. / FSmax ) / 2. );
                }
                if ( fabs ( FSN-FS ) <1.e-3 && conv==true ) {
                        anal.AcceptSolution();
                        //conv=true;
                }

                counterout++;
                if ( ( FSmax - FSmin ) / FS < tol  && conv==true ) {
                        anal.AcceptSolution();
                        //conv=true;
                }
        }  while ( ( ( FSmax - FSmin ) / FS > tol || conv==false ) && counterout<maxcout );

        auto t3 = chrono::high_resolution_clock::now();
        auto timeinmili = chrono::duration_cast<chrono::seconds> ( t3 - t0 );
        std::cout << "final safety factor "<< FS << " total time in ShearRed = "<< timeinmili.count() << " s "<<std::endl;
        return ( FSmax + FSmin ) /2;
}
//Tranfere a solucao nodal da malha cmesh para os pontos de integracao da malha fCompMesh. Este metodo e usado para transferir a solucao dos
// //campos estocasticos. O metodo findelement e caro, e custa muito ao monte carlo.
void SlopeAnalysis::TransferFieldsSolutionFrom ( int isol )
{

        TPZMatWithMem<TPZElastoPlasticMem> *pMatWithMem2 = dynamic_cast<TPZMatWithMem<TPZElastoPlasticMem> *> ( fCompMesh->MaterialVec() [1] );
        if ( pMatWithMem2 ) {
                pMatWithMem2->SetUpdateMem ( true );
        } else {
                DebugStop();
        }

        //malha geometrica field
        //const TPZGeoMesh * gmeshfield = fAnalysisField->Mesh()->Reference();

        //campos estocasticos
        TPZVec<TPZFMatrix<REAL>> fields=GetFields();

        int nfields = fields.size();

        //fields[0].Print("coes");

        for ( int imesh=0; imesh<nfields; imesh++ ) {
                // TPZCompMesh * mesh = new TPZCompMesh(*fCompMeshField->Clone());
                //vecfieldmesh[ifield] =fCompMeshField;
                //vecfieldmesh[ifield]->LoadSolution(fields[ifield]);

                fCompMeshField->LoadSolution ( fields[imesh] );

                //num elementos malha elastoplastica
                int nels =  fCompMesh->NElements();
                for ( int iel=0; iel<nels; iel++ ) {

                        TPZCompEl *cel = fCompMesh->ElementVec() [iel];
                        TPZInterpolationSpace *intel = dynamic_cast<TPZInterpolationSpace *> ( cel );
                        if ( !cel || !intel || dynamic_cast<TPZMatWithMem<TPZElastoPlasticMem> *> ( intel->Material() ) != pMatWithMem2 || ( intel->Material()->Id() <1 ) ) {
                                continue;
                        }


                        const TPZIntPoints &intpoints = intel->GetIntegrationRule();
                        int nint = intpoints.NPoints();
                        TPZManVector<REAL,3> point ( 2,0. );
                        // TPZVec<REAL> point ( 3,0. );


                        TPZMaterialDataT<REAL> data;
                        intel->InitMaterialData ( data );
                        data.fNeedsSol = true;

                        for ( long ip =0; ip<nint; ip++ ) {
                                REAL weight;
                                intpoints.Point ( ip, point, weight );
                                data.intLocPtIndex = ip;
                                intel->ComputeRequiredData ( data, point );

                                int indexplastic =data.intGlobPtIndex;
                                TPZElastoPlasticMem &mem = pMatWithMem2->MemItem ( indexplastic );
                                mem.m_elastoplastic_state.fmatprop.Resize ( 3 );
                                mem.m_elastoplastic_state.fmatpropinit.Resize ( 3 );


                                long elementid1 = 0;
                                TPZManVector<REAL,3> qsi ( 2,0. );
                                //TPZVec<REAL> qsi(3,0.);
                                int targetdim=2;
                                TPZCompEl*celfield;
                                TPZGeoEl *gelfield;

                                gelfield = fCompMeshField->Reference()->FindElement ( data.x, qsi, elementid1,targetdim );
                                celfield = gelfield->Reference();
                                if ( !celfield ) {
                                        std::cout << "Elemento computacional nao encontrado."<<endl;
                                        DebugStop();
                                        //continue;
                                }

                                TPZInterpolationSpace *intelfield = dynamic_cast<TPZInterpolationSpace *> ( celfield );



                                //cout <<"aaaaa"<<endl;
                                TPZMaterialDataT<REAL> datafield;
                                datafield.fNeedsSol = true;
                                //cout <<"bbbbb"<<endl;
                                intelfield->InitMaterialData ( datafield );

                                datafield.fNeedsSol = true;
                                intelfield->ComputeRequiredData ( datafield, qsi );
                                //cout <<"ccccc"<<endl;
                                REAL datarealvalue=datafield.sol[isol][0];

                                mem.m_elastoplastic_state.fmatpropinit[imesh] = datarealvalue;
                                mem.m_elastoplastic_state.fmatprop[imesh] = datarealvalue;

                        }
                        pMatWithMem2->SetUpdateMem ( false );

                }

        }
        //fCompMeshField->Solution().Redim(0,0);
}
void SlopeAnalysis::IntegrateFieldOverARegion(int imc)
{
        string saida = "post/regionmean";
        auto var=to_string ( imc );
        saida+=var;
        saida+=".dat";
        ofstream out ( saida );

        //campos estocasticos
        TPZVec<TPZFMatrix<REAL>> fields=GetFields();

        int nfields = fields.size();
        TPZVec<REAL> integrationmeanvalues(nfields);


        for ( int imesh=0; imesh<nfields; imesh++ ) {

                TPZFMatrix<REAL> solu(fields[imesh].Rows());
                for(int inodalsol=0;inodalsol< fields[imesh].Rows();inodalsol++)solu(inodalsol,0)=fields[imesh](inodalsol,imc);
                fCompMeshField->LoadSolution ( solu );

                //num elementos malha elastoplastica
                int nels =  fCompMeshField->NElements();
                REAL val=0.;
                int cout=0;
                for ( int iel=0; iel<nels; iel++ ) {

                        TPZCompEl *cel = fCompMeshField->ElementVec() [iel];
                        TPZInterpolationSpace *intel = dynamic_cast<TPZInterpolationSpace *> ( cel );
                        if (!cel || !intel) {
                                continue;
                        }

                        TPZManVector<REAL,3> point ( 2,0. );

                        TPZMaterialDataT<REAL> data;
                        intel->InitMaterialData ( data );
                        data.fNeedsSol = true;
                        intel->ComputeRequiredData(data,point);
                        REAL x=data.x[0];
                        REAL y=data.x[1] ;
                        REAL residual = (x-40.)*(x-40.)+(y-40.)*(y-40.)-20.*20.;
                        if(residual<0)
                        {
                            TPZVec<REAL> integral = cel->IntegrateSolution(0);
                            val+=integral[0];
                            cout++;
                        }

            }
            integrationmeanvalues[imesh]=val/cout;
        }
        out <<  integrationmeanvalues[0] <<endl;
        out << integrationmeanvalues[1] <<endl;
}

void SlopeAnalysis::ShearReductionIntegrationPoints ( REAL FS )
{
        TPZMatWithMem<TPZElastoPlasticMem> *pMatWithMem2 = dynamic_cast<TPZMatWithMem<TPZElastoPlasticMem> *> ( fCompMesh->MaterialVec() [1] );
        if ( pMatWithMem2 ) {
                pMatWithMem2->SetUpdateMem ( true );
        } else {
                DebugStop();
        }


        int nels =  fCompMesh->NElements();

        for ( int iel=0; iel<nels; iel++ ) {

                TPZCompEl *cel = fCompMesh->ElementVec() [iel];
                TPZInterpolationSpace *intel = dynamic_cast<TPZInterpolationSpace *> ( cel );
                if ( !cel || !intel || dynamic_cast<TPZMatWithMem<TPZElastoPlasticMem> *> ( intel->Material() ) != pMatWithMem2 || ( intel->Material()->Id() <1 ) ) {
                        continue;
                }


                const TPZIntPoints &intpoints = intel->GetIntegrationRule();
                int nint = intpoints.NPoints();
                TPZManVector<REAL,3> point ( 2,0. );


                TPZMaterialDataT<REAL> data;
                intel->InitMaterialData ( data );
                data.fNeedsSol = true;

                for ( long ip =0; ip<nint; ip++ ) {
                        REAL weight;
                        intpoints.Point ( ip, point, weight );
                        data.intLocPtIndex = ip;
                        intel->ComputeRequiredData ( data, point );

                        int indexplastic =data.intGlobPtIndex;
                        TPZElastoPlasticMem &mem = pMatWithMem2->MemItem ( indexplastic );
                        //mem.m_elastoplastic_state.fmatprop.Resize ( 3 );
                        if ( !mem.m_elastoplastic_state.fmatprop.size() ) {
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

TPZElastoPlasticAnalysis   SlopeAnalysis::SetSlopeAnalysis ( int type,int numthreads )
{
        TPZElastoPlasticAnalysis anal(fCompMesh,cout);

        switch ( type ) {
        case 0: {
                //cout << "Solver called with TPZStepSolver\n";
                TPZSkylineStructMatrix<STATE> matskl ( fCompMesh );
                matskl.SetNumThreads ( numthreads );
                anal.SetStructuralMatrix ( matskl );
                TPZStepSolver<STATE> step;
                step.SetDirect ( ELDLt );
                anal.SetSolver ( step );
                break;
        }
        case 1: {
                //cout << "Solver called with TPZPardisoSolver\n";
                TPZSSpStructMatrix<STATE> SSpStructMatrix ( fCompMesh );
                SSpStructMatrix.SetNumThreads ( numthreads );
                anal.SetStructuralMatrix ( SSpStructMatrix );
                TPZPardisoSolver<REAL> *pardiso = new TPZPardisoSolver<REAL>;
                anal.SetSolver ( *pardiso );
                break;
        }
        default: {
                cout << "Solver was not initialized properly\n";
                DebugStop();
        }
        }
        return anal;
}
void SlopeAnalysis::ManageFieldCretion()
{
        int nfields = fMeanvec.size();
        fFields.resize ( nfields );
        if ( !nfields ) DebugStop();
        for ( int ifield=0; ifield<nfields; ifield++ ) {
                fFields[ifield] = GenerateRandomField ( fMeanvec[ifield],fCovvec[ifield],fSamples );
        }

}
TPZFMatrix<REAL>  SlopeAnalysis::GenerateRandomField ( REAL mean, REAL cov,int samples )
{
        TPZFMatrix<REAL>  PHIt,PHI=fSolutionValVec;

        //std::cout <<" mean = "<< mean  << " cov = " << cov<<std::endl;
        if ( PHI.Rows() <=1 ) {
                std::cout <<" no solution"<< std::endl;
                DebugStop();
        }


        int M = PHI.Cols();

        TPZFMatrix<REAL> THETA ( M,samples );

        std::normal_distribution<REAL> distribution ( 0., 1. );
        for ( int n = 0; n < samples; n++ ) {
                for ( int iexp = 0; iexp < M; iexp++ ) {
                        std::random_device rd{};
                        std::mt19937 generator{ rd() };
                        REAL xic = distribution ( generator );
                        THETA ( iexp,n ) = xic;
                }
        }

        TPZFMatrix<REAL> hhat;
        PHI.Multiply ( THETA, hhat );

        REAL sdev = cov* mean;
        REAL xi = sqrt ( log ( 1 + pow ( ( sdev / mean ),2 ) ) );
        REAL lambda = log ( mean ) - xi * xi / 2.;

        for ( int i = 0; i < hhat.Rows(); i++ ) {
                for ( int j = 0; j < hhat.Cols(); j++ ) {
                        hhat ( i,j ) = exp ( lambda + xi * hhat ( i,j ) );
                }
        }

        //PrintField(hhat,filename);
        //LoadSolution(hhat);
        return hhat;
}

TPZCompMesh * SlopeAnalysis::CreateCMesh ( TPZGeoMesh *gmesh, int pOrder, REAL coes,REAL atrito )
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
        elasticresponse.SetEngineeringData ( E,nu );

        // Mohr Coulomb data
        REAL mc_cohesion    = coes;//kPa
        REAL mc_phi         = atrito;
        REAL mc_psi         = mc_phi;

        //elasticresponse.Print(std::cout);

        plasticmorh mohrcoulombplasticstep;

        mohrcoulombplasticstep.fYC.SetUp ( mc_phi, mc_psi, mc_cohesion, elasticresponse );

        mohrcoulombplasticstep.fER = elasticresponse;

        //mohrcoulombplasticstep.Print(std::cout);

        int PlaneStrain = 1;

        plasticmat * material = new plasticmat ( matid,PlaneStrain );

        material->SetPlasticityModel ( mohrcoulombplasticstep );

        material->SetId ( matid );

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
void SlopeAnalysis::DivideElementsAbove(REAL refineaboveval, std::set<long> &elindices)
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

void SlopeAnalysis::PRefineElementsAbove(REAL refineaboveval, int porder, std::set<long> &elindices)
{

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

void SlopeAnalysis::ComputeElementDeformation()
{
    long nelem = fCompMesh->NElements();
    fPlasticDeformSqJ2.resize(nelem);
    fPlasticDeformSqJ2.Fill(0.);
    fCompMesh->ElementSolution().Redim(nelem, 1);
    TPZMatWithMem<TPZElastoPlasticMem> *pMatWithMem2 = dynamic_cast<TPZMatWithMem<TPZElastoPlasticMem> *> (fCompMesh->MaterialVec()[1]);
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
                //pMatWithMem2->GetPlasticity().Phi(total,phi);
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

void SlopeAnalysis::PostPlasticity(std::string vtkd)
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

void  SlopeAnalysis::CreatePostProcessingMesh (TPZPostProcAnalysis * PostProcess )
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


void SlopeAnalysis::PostProcessVariables ( TPZStack<std::string> &scalNames, TPZStack<std::string> &vecNames )
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


TPZGeoMesh * SlopeAnalysis::TriGMesh ( int ref )
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
                                new TPZGeoElRefPattern< pzgeom::TPZGeoLinear> ( iel, TopoLine, -1, *gmesh );//bottom
                        } else if ( ( fabs ( ( x0-L ) ) <tol && fabs ( ( xf-L ) ) <tol ) ) {
                                new TPZGeoElRefPattern< pzgeom::TPZGeoLinear> ( iel, TopoLine, -2, *gmesh );//rigth
                        } else if ( ( fabs ( ( y0-h1 ) ) <tol && fabs ( ( yf-h1 ) ) <tol ) ) {
                                new TPZGeoElRefPattern< pzgeom::TPZGeoLinear> ( iel, TopoLine, -3, *gmesh );//toprigth
                        } else if ( ( fabs ( ( y0- ( h1+h2 ) ) ) <tol && fabs ( ( yf- ( h1+h2 ) ) ) <tol ) ) {
                                new TPZGeoElRefPattern< pzgeom::TPZGeoLinear> ( iel, TopoLine, -4, *gmesh );//topleft
                        } else if ( ( fabs ( ( x0-0 ) ) <tol && fabs ( ( xf-0 ) ) <tol ) ) {
                                new TPZGeoElRefPattern< pzgeom::TPZGeoLinear> ( iel, TopoLine, -5, *gmesh );//left
                        } else if ( ( fabs ( ( xf-x0 ) ) >tol && fabs ( ( yf-y0 ) ) >tol ) ) {
                                new TPZGeoElRefPattern< pzgeom::TPZGeoLinear> ( iel, TopoLine, -6, *gmesh );//ramp
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


void SlopeAnalysis::Write(TPZStream &buf, int withclassid) const
{
//     buf.Write(&fCohesion);
//     buf.Write(&fAtrito);
//     buf.Write(&fGammaW);
//     buf.Write(&fGammaS);
//     buf.Write(fMeanvec);
//     buf.Write(fCovvec);
//     buf.Write(fSamples);
//     buf.Write(fPlasticDeformSqJ2);
//     //fCompMeshField->Reference()->Write(buf, withclassid);
//     //fCompMeshField->Write(buf, withclassid);
//     fGMesh->Write(buf, 0);
//     fCompMesh->Write(buf, 0);
    fSolutionValVec.Write(buf,withclassid);
    fFields[0].Write(buf,withclassid);
    fFields[1].Write(buf,withclassid);
}

void SlopeAnalysis::Read(TPZStream &buf, void *context)
{
//     cout << "a"<<endl;
//     buf.Read(&fCohesion);
//     buf.Read(&fAtrito);
//     buf.Read(&fGammaW);
//     buf.Read(&fGammaS);
//     cout << "b"<<endl;
//     buf.Read(fMeanvec);
//     buf.Read(fCovvec);
//     buf.Read(&fSamples);
//     buf.Read(fPlasticDeformSqJ2);
//     cout << "c"<<endl;
//     //fCompMeshField->Reference()->Read(buf, context);
//     //fCompMeshField->Read(buf, context);
//     //fCompMesh->Reference()->Read(buf, context);
//     fGMesh->Read(buf, 0);
//     cout << "d"<<endl;
//     fCompMesh->Read(buf, fGMesh);
//     cout << "e"<<endl;
    fSolutionValVec.Read(buf,context);
    fFields.resize(2);
    fFields[0].Read(buf,context);
    fFields[1].Read(buf,context);
    cout << "f"<<endl;
}
  int SlopeAnalysis::ClassId() const
{
  return Hash("SlopeAnalysis");
}
