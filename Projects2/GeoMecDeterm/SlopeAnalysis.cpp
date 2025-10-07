// SPDX-FileCopyrightText: 2024 <copyright holder> <email>
// SPDX-License-Identifier: Apache-2.0

#include "SlopeAnalysis.h"

SlopeAnalysis::SlopeAnalysis()
        : fCompMesh ( 0 ), fGMesh ( 0 ),  fNumThreads ( 0 ),
          fSolver ( 0 )
{
        // Construtor padrão
}

SlopeAnalysis::SlopeAnalysis ( const SlopeAnalysis& other )
:  fCompMesh ( other.fCompMesh ), fGMesh ( other.fGMesh  ),  fNumThreads ( other.fNumThreads  ),
fSolver ( other.fSolver)
{


}
SlopeAnalysis::SlopeAnalysis ( TPZGeoMesh * gmesh, TPZCompMesh * cmesh,int nThreads,int solver)
:  fCompMesh ( cmesh ), fGMesh ( gmesh ),  fNumThreads ( nThreads),fSolver ( solver)
{

}

SlopeAnalysis::~SlopeAnalysis()
{
        // Destrutor, limpando a memória alocada
        delete fCompMesh;
        delete fGMesh;

}


REAL SlopeAnalysis::IterativeProcessArcLength ( REAL tol,int numiter,REAL tol2,int numiter2,REAL l,REAL lambda0,bool &converge )
{


        TPZElastoPlasticAnalysis anal =  SetSlopeAnalysis ( );

        std::vector<double> fslist;

        REAL lambda=lambda0;

        TPZFMatrix<REAL> u,dws,dwb,du,rhs,rhs1,rhs2,uold;

        u=anal.Solution();
        u.Zero();

        LoadingRamp ( 1. );
        anal.AssembleResidual();
        TPZFMatrix<REAL> rhstotal=anal.Rhs();

        LoadingRamp ( 0. );
        anal.AssembleResidual();
        TPZFMatrix<REAL> rhsint=anal.Rhs();
        TPZFMatrix<REAL> FBODY=rhstotal-rhsint;
        REAL normint=Norm(rhstotal);
        REAL diff=1000,lambdan;
        int counterout=0;
        REAL ndesi=7;
        REAL fac=1;
        REAL l0=l;
        do {
                cout << "\n load step  = " << counterout+1 << " load factor  = " << lambda << " diff = " << diff << " l = " << l<< " fac = "<< fac <<  endl;
                int counter=0;
                REAL normrhs=10.;
                REAL normrhsold=10.;

                REAL dlamb=0.;
                uold=u;
                //anal.Solution().Zero();

                lambdan=lambda;
                //lambda=lambda0 ;


                u.Zero();
                anal.LoadSolution ( u );
                do { //while( counter<numiter2 && normdu>tol2 );


                        //K dws = -r
                        LoadingRamp ( lambda );
                        anal.Assemble();
                        anal.Solve();
                        dws=anal.Solution();
                        rhs=anal.Rhs();

                        anal.Rhs() =FBODY;
                        anal.Solve();
                        dwb=anal.Solution();


                        normrhsold=normrhs;
                        normrhs=Norm ( rhs )/normint;


                        cout << "i = " << counter  << " ||rhs|| =  "<< normrhs <<" lambda = "<< lambda << " dlamb = "<< dlamb <<" l = "<< l <<endl;

                        TPZVec<REAL> lambvec;
                        if ( counter == 0 ) {
                                dlamb =computelamda0 ( dwb, u, l );
                        } else {

                                if ( false ) {

                                        dlamb = computelamda ( dwb, dws, u, l );
                                } else {
                                        TPZVec<REAL> lambvec = computelamdacris ( dwb, dws, u, l );

                                        LoadingRamp ( lambvec[0]+lambda );
                                        anal.AssembleResidual();
                                        rhs1=anal.Rhs();

                                        LoadingRamp ( lambvec[1]+lambda );
                                        anal.AssembleResidual();
                                        rhs2=anal.Rhs();

                                        REAL res1 = Norm ( rhs1 );
                                        REAL res2 = Norm ( rhs2 );

                                        if ( res1 <res2 ) {
                                                dlamb = lambvec[0];

                                        } else {
                                                dlamb = lambvec[1];
                                        }
                                }

                        }



                        if(normrhs>normrhsold)
                        {
                                //lambda-=dlamb;
                                //normrhs=normrhsold;
                                //dlamb=0.;
                                //u=uold;
                                //lambda=lambda0 ;
                             //   l*=0.5;
                                //u=uold;
                                //continue;
                        }
                        //else{

                                lambda += dlamb;
                                du=dws+dlamb*dwb;
                                //uold=u;
                                u+=du;
                                anal.LoadSolution ( u );
                       // }


                        counter++;

                } while ( counter<numiter2 && normrhs>tol2 );

                counterout++;

                anal.AcceptSolution();
                diff=fabs ( lambda-lambdan );


                fac=ndesi  / ( counter+1 );

                //l=fac*l;
                if(l<0.01)
                {
                        l=0.01;
                }
                if(l>5)
                {
                        l=5.;
                }

                cout << "fac = " << fac << "(counter+1)  = "<< ( counter+1 ) << "lold = "<<l << " diff = "<< diff<<endl;

        } while ( counterout<numiter && diff>tol );

        if(diff>tol)
        {
                converge=false;
                cout << "Convergence in arc length faile whith diff = "<< diff << endl;
        }
        else
        {
                cout << "arc length converged whith diff = "<< diff << endl;
                converge=true;
        }

        anal.AcceptSolution();
        return lambda;
}

REAL SlopeAnalysis::ArcLength(bool &conv)
{
        REAL tol=0.01;
        int numiter=200;
        REAL tol2=0.001;
        int numiter2=60;
        REAL l=0.5;
        REAL lambda0=1;

        REAL FS  = IterativeProcessArcLength ( tol,numiter,tol2,numiter2,l,lambda0,conv );

        return FS;
}

REAL SlopeAnalysis::SolveDeterministic ( bool IsSRM ,REAL coes,REAL atrito)
{
        REAL FSOLD,FS;
        int neqold;
        int neq=fCompMesh->NEquations();
        cout << "NUMBER OF EQUATIONS  = " << neq << endl;
        InitializeMemory(coes,atrito);
        bool conv;

        if ( IsSRM==true ) {
                //FS = ShearRed ( 20,0.5,0.01 );
                FS =ShearRedNoIntegrationPoints( 20,0.5,0.01 );
        } else {
                //FS = GravityIncrease() ;

                FS  =ArcLength(conv);
                if(conv==false)
                {
                   FS = GravityIncrease() ;
                }
        }

        cout << "Refining.."<<endl;
        std::set<long> elindices,elindices2;
        for ( int iref=1; iref<=1; iref++ ) {
                cout << "computing deformation..."  << endl;
                ComputeElementDeformation();
                cout << "p refining..."  << endl;
                PRefineElementsAbove ( 0.01, fCompMesh->GetDefaultOrder()+iref,elindices2 );
                cout << "h refining..."  << endl;
                DivideElementsAbove ( 0.01,elindices );
                neqold=neq;
                neq=fCompMesh->NEquations();
                cout << "initializing memory..."  << endl;
                InitializeMemory(coes,atrito);
                cout << "# of equations  = " <<neq << " fabs(FS-FSOLD)  "  << fabs ( FS-FSOLD )  << endl;
                FSOLD=FS;
                if ( IsSRM==true ) {
                       // FS = ShearRed ( 20,FSOLD,0.01 );
                        FS =ShearRedNoIntegrationPoints( 20,FSOLD,0.01 );
                } else {

                //FS = ShearRed ( 20,FSOLD,0.01 );
                                        FS  =ArcLength(conv);
                if(conv==false)
                {
                   FS = GravityIncrease() ;
                }
                }

                if ( fabs ( FS-FSOLD ) <0.01 ) {
                        //cout << " FS-FSOLD = "<< fabs ( FS-FSOLD ) <<endl;
                        //break;
                } else if ( neq>10000||neqold==neq ) {
                        //cout << " neq>10000 = "<< neq <<" neqold = "<< neqold <<endl;
                        //break;
                }
        }

        InitializeMemory(coes,atrito);
        string meshref = "refinidemesh-grid";
        meshref+=".vtk";
        std::ofstream files ( meshref );
        TPZVTKGeoMesh::PrintGMeshVTK ( fCompMesh->Reference(),files,true );

        return FS;
}


void SlopeAnalysis::ApplyGravityLoad ( TPZManVector<REAL, 3> bodyforce )
{
        plasticmat * body= dynamic_cast<plasticmat *> ( fCompMesh->FindMaterial ( 1 ) );
        // body->SetBodyForce ( bodyforce );
        //REIMPLEMENTAR
        DebugStop();
        //body->SetLoadFactor ( factor );
}

void SlopeAnalysis::LoadingRamp ( REAL factor )
{
        plasticmat * body= dynamic_cast<plasticmat *> ( fCompMesh->FindMaterial ( 1 ) );
        TPZManVector<REAL, 3> force ( 3,0. );
        //REIMPLEMENTAR
        DebugStop();
        //body->SetLoadFactor ( factor );


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

                TPZElastoPlasticAnalysis anal =  SetSlopeAnalysis ( );
                fCompMesh->Solution().Zero();
                REAL norm = 1000.;
                REAL tol2 = 1.e-3;
                int NumIter = 100;
                bool linesearch = true;
                bool checkconv = false;
                int iters;

                ShearReductionIntegrationPoints ( FS );

                auto t1 = chrono::high_resolution_clock::now();
                conv  =anal.IterativeProcess ( cout,tol2, NumIter,  linesearch,  checkconv,iters );

//                 int numit1=20,numit2=10;
//                 REAL tolfs=1.e-2,tolrhs=1.e-3,l=0.5;
//                 FS = IterativeProcessArcLength ( tolfs,numit1,tolrhs,numit2,l,FS,conv );


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

REAL SlopeAnalysis::ShearRedNoIntegrationPoints ( int maxcout,REAL FS0,REAL fstol )
{
        LoadingRamp ( 1. );

        plasticmat * body= dynamic_cast<plasticmat *> ( fCompMesh->FindMaterial ( 1 ) );
        REAL FS=FS0,FSmax=10.,FSmin=0.,tol=fstol;
        int counterout = 0;
        bool conv=false;
        auto t0 =chrono::high_resolution_clock::now();
        REAL FSN=1000;


        do {

                TPZElastoPlasticAnalysis anal =  SetSlopeAnalysis ( );

                body->GetPlasticModel().SetStrengthReductionFactor(FS);

                fCompMesh->Solution().Zero();
                REAL norm = 1000.;
                REAL tol2 = 1.e-3;
                int NumIter = 100;
                bool linesearch = true;
                bool checkconv = false;
                int iters;

                auto t1 = chrono::high_resolution_clock::now();
                conv  =anal.IterativeProcess ( cout,tol2, NumIter,  linesearch,  checkconv,iters );

//                 int numit1=20,numit2=10;
//                 REAL tolfs=1.e-2,tolrhs=1.e-3,l=0.5;
//                 FS = IterativeProcessArcLength ( tolfs,numit1,tolrhs,numit2,l,FS,conv );


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

        body->GetPlasticModel().SetStrengthReductionFactor(1.);
        auto t3 = chrono::high_resolution_clock::now();
        auto timeinmili = chrono::duration_cast<chrono::seconds> ( t3 - t0 );
        std::cout << "final safety factor "<< FS << " total time in ShearRed = "<< timeinmili.count() << " s "<<std::endl;
        return ( FSmax + FSmin ) /2;
}

REAL SlopeAnalysis::GravityIncrease ( )
{

        REAL FS=0.1,FSmax=1000.,FSmin=0.,tol=0.01;
        int neq = fCompMesh->NEquations();
        int maxcount=100;
        TPZFMatrix<REAL> displace ( neq,1 ),displace0 ( neq,1 );

        int counterout = 0;
        REAL factor =1.;
        LoadingRamp ( factor );
        REAL norm = 1000.;
        REAL tol2 = 1.e-2;
        int NumIter = 100;
        bool linesearch = true;
        bool checkconv = false;

        TPZElastoPlasticAnalysis anal =  SetSlopeAnalysis ( );

        anal.AssembleResidual();

         TPZFMatrix<REAL> fbody=anal.Rhs();
         REAL normfbofy=Norm(fbody);

         cout << " norm fbody = "<< normfbofy << endl;

        do {

                std::cout << "FS = " << FS  <<" | Load step = " << counterout << " | Rhs norm = " << norm  << std::endl;
                LoadingRamp ( FS );
                //SetSuportPressure(cmesh,FS);

                TPZElastoPlasticAnalysis anal =  SetSlopeAnalysis ( );
                chrono::steady_clock sc;
                auto start = sc.now();
                int iters;
                bool conv  =anal.IterativeProcess ( cout,tol2, NumIter,  linesearch,  checkconv,iters );
                //bool conv =anal->IterativeProcess(cout, tol2, NumIter,linesearch,checkconv);
                auto end = sc.now();
                auto time_span = static_cast<chrono::duration<double>> ( end - start );
                //cout << "| total time in iterative process =  " << time_span.count()<< std::endl;
                //anal->IterativeProcess ( outnewton, tol2, NumIter);

                norm = Norm ( anal.Rhs() )/normfbofy;

                if ( conv==false ) {
                        fCompMesh->LoadSolution ( displace0 );
                        //cmesh->Solution().Zero();
                        FSmax = FS;
                        FS = ( FSmin + FSmax ) / 2.;

                } else {
                        // uy+=findnodalsol(cmesh);
                        displace0 = anal.Solution();
                        FSmin = FS;
                        anal.AcceptSolution();
                        FS = 1. / ( ( 1. / FSmin + 1. / FSmax ) / 2. );
                }
                // cout << "|asdadadasd =  " << std::endl;
                counterout++;

        }  while ( ( ( FSmax - FSmin ) / FS > tol && counterout<maxcount ) );


        return FS;
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

TPZElastoPlasticAnalysis   SlopeAnalysis::SetSlopeAnalysis ( )
{
        TPZElastoPlasticAnalysis anal ( fCompMesh,cout );

        switch ( fSolver ) {
        case 0: {
                //cout << "Solver called with TPZStepSolver\n";
                TPZSkylineStructMatrix<STATE> matskl ( fCompMesh );
                matskl.SetNumThreads ( fNumThreads );
                anal.SetStructuralMatrix ( matskl );
                TPZStepSolver<STATE> step;
                step.SetDirect ( ELDLt );
                anal.SetSolver ( step );
                break;
        }
        default: {
                cout << "Solver was not initialized properly\n";
                DebugStop();
        }
        }
        return anal;
}



void SlopeAnalysis::DivideElementsAbove ( REAL refineaboveval, std::set<long> &elindices )
{
        //int porder =fPorder+3;
        //fGmesh->ResetReference();
        //fCompMesh->LoadReferences();
        TPZManVector<REAL,3> findel ( 3,0. ),qsi ( 2,0. );


        long nelem = fCompMesh->NElements();
        for ( long el=0; el<nelem; el++ ) {
                TPZCompEl *cel = fCompMesh->ElementVec() [el];
                if ( !cel ) {
                        continue;
                }

                TPZInterpolationSpace *intel = dynamic_cast<TPZInterpolationSpace *> ( cel );
                if ( !intel ) {
                        DebugStop();
                }
                TPZMatWithMem<TPZElastoPlasticMem> *pMatWithMem2 = dynamic_cast<TPZMatWithMem<TPZElastoPlasticMem> *> ( cel->Material() );
                if ( !pMatWithMem2 ) {
                        continue;
                }
                const TPZFMatrix<STATE> &elsol = fCompMesh->ElementSolution();
                if ( elsol.Get ( el,0 ) <=refineaboveval ) {
                        continue;
                }
                //intel->PRefine(3);
                int porder = intel->GetPreferredOrder();
                TPZStack<long> subels;
                long index = cel->Index();


                intel->Divide ( index, subels,0 );
                for ( int is=0; is<subels.size(); is++ ) {
                        elindices.insert ( subels[is] );
                        TPZCompEl *subcel = fCompMesh->ElementVec() [subels[is]];

                        TPZInterpolationSpace *subintel = dynamic_cast<TPZInterpolationSpace *> ( subcel );
                        if ( !subintel ) {
                                DebugStop();
                        }
                        subintel->SetPreferredOrder ( porder );
                }
        }
        // divide elements with more than one level difference
        bool changed = true;
        while ( changed ) {
                changed = false;
                std::set<long> eltodivide;
                long nelem = fCompMesh->NElements();
                for ( long el=0; el<nelem; el++ ) {
                        TPZCompEl *cel = fCompMesh->ElementVec() [el];
                        if ( !cel ) {
                                continue;
                        }
                        TPZInterpolationSpace *intel = dynamic_cast<TPZInterpolationSpace *> ( cel );
                        if ( !intel ) {
                                DebugStop();
                        }
                        TPZGeoEl *gel = cel->Reference();
                        if ( !gel ) {
                                DebugStop();
                        }
                        int ns = gel->NSides();
                        for ( int is=0; is<ns; is++ ) {
                                TPZGeoElSide gelside ( gel, is );
                                if ( gelside.Dimension() != 1 ) {
                                        continue;
                                }
                                TPZCompElSide big = gelside.LowerLevelCompElementList2 ( 1 );
                                if ( !big ) {
                                        continue;
                                }
                                TPZGeoElSide geobig ( big.Reference() );
                                // boundary elements will be refined by AdjustBoundaryElements
                                if ( geobig.Element()->Dimension() != 2 ) {
                                        continue;
                                }
                                if ( gel->Level()-geobig.Element()->Level() > 1 ) {
                                        eltodivide.insert ( big.Element()->Index() );
                                }
                        }
                }
                std::set<long>::iterator it;
                for ( it = eltodivide.begin(); it != eltodivide.end(); it++ ) {
                        changed = true;
                        long el = *it;
                        TPZCompEl *cel = fCompMesh->ElementVec() [el];
                        if ( !cel ) {
                                continue;
                        }
                        TPZInterpolationSpace *intel = dynamic_cast<TPZInterpolationSpace *> ( cel );
                        if ( !intel ) {
                                DebugStop();
                        }
                        TPZMatWithMem<TPZElastoPlasticMem> *pMatWithMem2 = dynamic_cast<TPZMatWithMem<TPZElastoPlasticMem> *> ( cel->Material() );
                        if ( !pMatWithMem2 ) {
                                continue;
                        }

                        int porder = intel->GetPreferredOrder();
                        TPZStack<long> subels;
                        long index = cel->Index();
                        intel->Divide ( index, subels,0 );
                        for ( int is=0; is<subels.size(); is++ ) {
                                elindices.insert ( subels[is] );
                                TPZCompEl *subcel = fCompMesh->ElementVec() [subels[is]];
                                TPZInterpolationSpace *subintel = dynamic_cast<TPZInterpolationSpace *> ( subcel );
                                if ( !subintel ) {
                                        DebugStop();
                                }
                                subintel->SetPreferredOrder ( porder );
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

void SlopeAnalysis::PRefineElementsAbove ( REAL refineaboveval, int porder, std::set<long> &elindices )
{

        fCompMesh->LoadReferences();
        long nelem = fCompMesh->NElements();
        for ( long el=0; el<nelem; el++ ) {
                TPZCompEl *cel = fCompMesh->ElementVec() [el];
                if ( !cel ) {
                        continue;
                }
                TPZInterpolationSpace *intel = dynamic_cast<TPZInterpolationSpace *> ( cel );
                if ( !intel ) {
                        DebugStop();
                }
                TPZMatWithMem<TPZElastoPlasticMem> *pMatWithMem2 = dynamic_cast<TPZMatWithMem<TPZElastoPlasticMem> *> ( cel->Material() );
                if ( !pMatWithMem2 ) {
                        continue;
                }

                const TPZFMatrix<STATE> &elsol = fCompMesh->ElementSolution();
                if ( elsol ( el,0 ) < refineaboveval ) {
                        continue;
                }
                //cout << "porder = " << porder << endl;
                TPZStack<long> subels;
                long index = cel->Index();
                elindices.insert ( index );
                intel->SetPreferredOrder ( porder );
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
        fPlasticDeformSqJ2.resize ( nelem );
        fPlasticDeformSqJ2.Fill ( 0. );
        fCompMesh->ElementSolution().Redim ( nelem, 1 );
        TPZMatWithMem<TPZElastoPlasticMem> *pMatWithMem2 = dynamic_cast<TPZMatWithMem<TPZElastoPlasticMem> *> ( fCompMesh->MaterialVec() [1] );
        if ( !pMatWithMem2 ) {
                fPlasticDeformSqJ2.Fill ( 0. );
        } else {
                for ( long el = 0; el<nelem; el++ ) {
                        TPZCompEl *cel = fCompMesh->ElementVec() [el];
                        fPlasticDeformSqJ2[el] = 0.;
                        if ( !cel ) {
                                continue;
                        }
                        TPZManVector<long> memindices;
                        cel->GetMemoryIndices ( memindices );
                        int numind = memindices.size();
                        REAL sqj2el = 0.00;
                        for ( int ind=0; ind<numind; ind++ ) {
                                int memoryindex = memindices[ind];
                                if ( memoryindex < 0 ) {
                                        continue;
                                }
                                TPZElastoPlasticMem &mem = pMatWithMem2->MemItem ( memindices[ind] );
                                TPZTensor<REAL> plastic =mem.m_elastoplastic_state.EpsP();

                                REAL J2 = plastic.J2();
                                REAL sqj2 = sqrt ( J2 );

                                sqj2el = max ( sqj2,sqj2el );


                        }
                        fPlasticDeformSqJ2[el] = sqj2el;
                }
        }
        fCompMesh->SetElementSolution ( 0, fPlasticDeformSqJ2 );
}


void SlopeAnalysis::InitializeMemory (REAL coesion, REAL atrito )
{

        TPZMatWithMem<TPZElastoPlasticMem> *pMatWithMem2 = dynamic_cast<TPZMatWithMem<TPZElastoPlasticMem> *> ( fCompMesh->MaterialVec() [1] );
        if ( pMatWithMem2 ) {
                pMatWithMem2->SetUpdateMem ( true );
        } else {
                DebugStop();
        }

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

                        mem.m_elastoplastic_state.fmatpropinit[0] = coesion;
                        mem.m_elastoplastic_state.fmatpropinit[1] = atrito;
                        //asscociativo
                        mem.m_elastoplastic_state.fmatpropinit[2] = atrito;

                        mem.m_elastoplastic_state.fmatprop[0] = coesion;
                        mem.m_elastoplastic_state.fmatprop[1] = atrito;
                        //asscociativo
                        mem.m_elastoplastic_state.fmatprop[2] = atrito;

                }

        }

        pMatWithMem2->SetUpdateMem ( false );

}


void SlopeAnalysis::PostPlasticity ( std::string vtkd )
{
        TPZPostProcAnalysis * postprocdeter = new TPZPostProcAnalysis();
        CreatePostProcessingMesh ( postprocdeter );


        TPZVec<int> PostProcMatIds ( 1,1 );

        TPZStack<std::string> PostProcVars, scalNames, vecNames;

        PostProcessVariables ( scalNames, vecNames );

        //string vtkd = "postprocessdeter.vtk";
        postprocdeter->DefineGraphMesh ( 2,scalNames,vecNames,vtkd );

        postprocdeter->PostProcess ( 0 );

        auto var=vtkd;
        std::ofstream files ( vtkd );
        TPZVTKGeoMesh::PrintGMeshVTK ( fCompMesh->Reference(),files,true );

        delete postprocdeter;
}

void  SlopeAnalysis::CreatePostProcessingMesh ( TPZPostProcAnalysis * PostProcess )
{
        if ( PostProcess->ReferenceCompMesh() != fCompMesh ) {

                PostProcess->SetCompMesh ( fCompMesh );

                TPZVec<int> PostProcMatIds ( 1,1 );
                TPZStack<std::string> PostProcVars, scalNames, vecNames;
                PostProcessVariables ( scalNames, vecNames );

                for ( int i=0; i<scalNames.size(); i++ ) {
                        PostProcVars.Push ( scalNames[i] );
                }
                for ( int i=0; i<vecNames.size(); i++ ) {
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


int SlopeAnalysis::ClassId() const
{
        return Hash ( "SlopeAnalysis" );
}
REAL  SlopeAnalysis::computelamda0 ( TPZFMatrix<REAL>& dwb,  TPZFMatrix<REAL>& fext, REAL& l )
{

        TPZFMatrix<REAL> dwt,solsig;
        fext.Transpose ( &dwt );
        dwt.Multiply ( dwb,solsig );
        REAL scal = solsig.Get ( 0,0 );

        REAL signum=0;
        //page 111, eq. 4.123 - Souza Neto //verificar sinal
        if ( scal<0 ) {
                signum=-1;
        } else {
                signum=1;
        }


        TPZFMatrix<REAL> dwbt, aparam;
        dwb.Transpose ( &dwbt );
        dwbt.Multiply ( dwb, aparam );

        return signum*l/sqrt ( aparam.Get ( 0,0 ) ) ;


}

REAL  SlopeAnalysis::computelamda ( TPZFMatrix<REAL>& dwb, TPZFMatrix<REAL>& dws, TPZFMatrix<REAL>& dw, REAL& l )
{



        int sz = dwb.Rows();
        REAL aa = 0.;

        aa = Dot ( dwb,dwb );

        REAL bb = 0.;

        TPZFMatrix<REAL> dwcopy;

        dwcopy = dw+dws;

        bb = Dot ( dwb,dwcopy );

        bb *= 2;
        REAL cc = 0.;

        cc= Dot ( dwcopy,dwcopy );

        cc -= l * l;
        REAL delta = bb * bb - 4. * aa * cc;
        REAL dlamb2;
        REAL dlamb1;


        //cout << "delta = " << delta << endl;
        //cout << "aa = " << aa << endl;
        //cout << "bb = " << bb << endl;
        //cout << "cc = " << cc << endl;
        if ( fabs ( aa ) >1.e-12 && delta>0 ) {
                dlamb2 = ( -bb + sqrt ( delta ) ) / ( 2. * aa ); //maior
                dlamb1= ( -bb - sqrt ( delta ) ) / ( 2. * aa ); //menor
                //return dlamb1;
                //cout << "dlamb1" <<dlamb1 << " dlamb2 = "<< dlamb2 << endl;
        } else {
                if ( bb != 0 ) {
                        //cout << "-cc/bb" <<-cc/bb << endl;
                        return -cc/bb;
                } else {
                        //cout << "(-bb ) / (2. * aa)" <<(-bb ) / (2. * aa)<< endl;
                        return ( -bb ) / ( 2. * aa );
                }
        }


        //page 111, eq. 4.118 - Souza Neto
        TPZFMatrix<REAL> temp1,temp1t,sol1,temp2,temp2t,sol2;
        temp1=dwb;
        temp1*=dlamb1;
        temp1+=dws;
        temp1+=dw;
        temp1.Transpose ( &temp1t );
        temp1t.Multiply ( dw,sol1 );

        temp2=dwb;
        temp2*=dlamb2;
        temp2+=dws;
        temp2+=dw;
        temp2.Transpose ( &temp2t );
        temp2t.Multiply ( dw,sol2 );

        if ( sol1.Get ( 0,0 ) >=sol2.Get ( 0,0 ) ) {
                //cout << "return 1 " << " sol1.Get ( 0,0 ) "<< sol1.Get ( 0,0 ) << " sol2.Get ( 0,0 ) "<< sol2.Get ( 0,0 ) <<endl;
                return dlamb1;
        } else {
                //cout << "return 2 " << endl;
                return dlamb2;
        }
}
TPZVec<REAL> SlopeAnalysis::computelamdacris ( TPZFMatrix<REAL>& dwb, TPZFMatrix<REAL>& dws, TPZFMatrix<REAL>& dw, REAL& l )
{
        // Tamanho da matriz
        int sz = dwb.Rows();

        // Cálculo de 'aa'
        REAL aa = Dot ( dwb, dwb );

        // Cálculo de 'bb'
        TPZFMatrix<REAL> dwcopy = dw + dws;
        REAL bb = 2.0 * Dot ( dwb, dwcopy );

        // Cálculo de 'cc'
        REAL cc = Dot ( dwcopy, dwcopy ) - l * l;

        // Delta da equação quadrática
        REAL delta = bb * bb - 4.0 * aa * cc;

        // Vetor de lambdas
        TPZVec<REAL> lambvec ( 2, 0 );

        // Tratamento do caso delta >= 0
        if ( fabs ( aa ) > 1.e-12 && delta >= 0 ) {
                REAL sqrtDelta = sqrt ( delta );
                REAL inv2a = 1.0 / ( 2.0 * aa ); // Evita cálculo redundante

                REAL dlamb1 = ( -bb - sqrtDelta ) * inv2a; // Menor raiz
                REAL dlamb2 = ( -bb + sqrtDelta ) * inv2a; // Maior raiz

                lambvec[0] = dlamb2;
                lambvec[1] = dlamb1;

        } else if ( fabs ( bb ) > 1.e-12 ) {
                // Caso especial onde aa é pequeno e bb != 0
                lambvec[0] = -cc / bb;
        } else {
                // Caso degenerado
                lambvec[0] = -bb / ( 2.0 * aa );
        }

        // Retorna os valores calculados
        return lambvec;
}
