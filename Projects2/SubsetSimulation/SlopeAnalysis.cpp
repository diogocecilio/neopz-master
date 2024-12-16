// SPDX-FileCopyrightText: 2024 <copyright holder> <email>
// SPDX-License-Identifier: Apache-2.0

#include "SlopeAnalysis.h"

SlopeAnalysis::SlopeAnalysis()
        : fCohesion ( 0 ), fAtrito ( 0 ), fGammaW ( 0 ), fGammaS ( 0 ),
          fNSamples ( 0 ),fCompMesh ( 0 ), fGMesh ( 0 ),  fNumThreads ( 0 ),
          fSolver ( 0 ), fFieldSamples (0),fMeanvec(0),fCovvec(0),fFieldSamplesSubSetFN()
{
        // Construtor padrão
}

SlopeAnalysis::SlopeAnalysis ( const SlopeAnalysis& other )
        : fCohesion ( other.fCohesion ), fAtrito ( other.fAtrito ), fGammaW ( other.fGammaW ), fGammaS ( other.fGammaS ),
          fNSamples ( other.fNSamples ),fCompMeshField ( other.fCompMeshField ), fSolutionValVec ( other.fSolutionValVec ),fFields ( other.fFields ),fRef0 ( other.fRef0 ),fPorder ( other.fPorder ),  fNumThreads ( other.fNumThreads ),
          fSolver ( other.fSolver ), fFieldSamples (other.fFieldSamples),fMeanvec(other.fMeanvec),fCovvec(other.fCovvec),fFieldSamplesSubSetFN(other.fFieldSamplesSubSetFN)
{
        fGMesh = TriGMesh ( fRef0 );
        //fGMesh = QuadGMesh ( fRef0 );
        fCompMesh = CreateCMesh ( fGMesh, fPorder, fCohesion, fAtrito );
        SetSlopeAnalysis ( );

}

SlopeAnalysis::SlopeAnalysis ( REAL gammaagua, REAL gammasolo, REAL coes, REAL atrito, int ref0, int porder,int therads,int solver )
        : fCohesion ( coes ), fAtrito ( atrito ), fGammaW ( gammaagua ), fGammaS ( gammasolo ),fRef0 ( ref0 ),fPorder ( porder ),fNumThreads ( therads ),
          fSolver ( solver )
{
        fGMesh = TriGMesh ( fRef0 );
        //fGMesh = QuadGMesh ( fRef0 );
        fCompMesh = CreateCMesh ( fGMesh, fPorder, fCohesion, fAtrito );
        SetSlopeAnalysis ( );
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

        do {
                cout << "\n load step  = " << counterout+1 << " load factor  = " << lambda << " diff = " << diff << " l = " << l<< " fac = "<< fac <<  endl;
                int counter=0;
                REAL normrhs=10.;
                REAL normdu=10.;
                REAL normrhsold=10.;

                REAL dlamb=0.;
                uold=u;
                //anal.Solution().Zero();

                lambdan=lambda;
                lambda=lambda0 ;

                //diff=1000;
                u.Zero();
                anal.LoadSolution ( u );
                do {
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

                lambda += dlamb;
                du=dws+dlamb*dwb;
                normdu=Norm(du);
                u+=du;
                anal.LoadSolution ( u );

                cout << "i = " << counter  << " ||rhs|| =  "<< normrhs <<" ||du||  = " <<normdu <<" lambda = "<< lambda << " dlamb = "<< dlamb <<" l = "<< l <<endl;

                counter++;

                } while ( counter<numiter2 && (normdu>tol2 ||normrhs>tol2*10));

                counterout++;
                anal.AcceptSolution();
                diff=fabs ( lambda-lambdan );
                fac=ndesi  / ( counter+1 );

        //} while ( counterout<numiter );
        } while ( counterout<numiter && diff>tol );

        //anal.AcceptSolution();
        if(counterout < numiter)
        {
                cout << "Converged"<<endl;
                converge=true;
        }else{
                cout << "Not Converged"<<endl;
                converge=false;
        }

        return lambda;
}



void SlopeAnalysis::FindRoot (bool &conv )
{


        TPZElastoPlasticAnalysis anal =  SetSlopeAnalysis ( );
         plasticmat * body= dynamic_cast<plasticmat *> ( fCompMesh->FindMaterial ( 1 ) );


  REAL normrhs=10000,normdu=10000,normrhsn=10000,normdun=10000,normrhs0;
    TPZFMatrix<STATE> x(anal.Solution()), dx(anal.Solution());
    x.Zero();
    dx.Zero();
    REAL tol = 1.e-3;
    int n_it = 100;
    anal.AssembleResidual();
    normrhs0=Norm(anal.Rhs());

    for (int i = 1; i <= n_it; i++) {
        anal.Assemble();
        anal.Solve();
        if ( false)
        {
            TPZFMatrix<STATE> nextSol;
            REAL LineSearchTol = 0.001 * Norm ( anal.Solution() );
            const int niter =2;
            anal.LineSearch ( x, anal.Solution(), nextSol, LineSearchTol, niter );
            x = nextSol;
            dx=anal.Solution();
        }else{
        dx = anal.Solution();
        x += dx;
        }
        //body->GetPlasticModel().SetStrengthReductionFactor(lambda);
        anal.LoadSolution(x);

        normdun=normdu;
        normdu=Norm(dx);

        normrhsn=normrhs;
        normrhs = Norm(anal.Rhs())/normrhs0;

//std::cout <<"iter = "<< i << " normrhs= " << normrhs<< " normrhsn= " << normrhsn<<" normdu= " << normdu<< " normdun= " << normdun<< std::endl;
        if (normrhs<tol) {

                std::cout <<"iter = "<< i << " normrhs= " << normrhs<< " normrhsn= " << normrhsn<<" normdu= " << normdu<< " normdun= " << normdun<< std::endl;
                //anal.AcceptSolution();

                conv=true;
                return;
        }else if(normrhsn<normrhs&&normdun<normdu){
            std::cout << "Fail to converge. Divergent method." << std::endl;
            conv= false;
            return;
        }
    }

    std::cout << " Not converged. Maximum number of iterations reached." << std::endl;
    conv= false;
    return;
}


REAL SlopeAnalysis::ArcLength(bool &conv)
{
        REAL tol=0.01;
        int numiter=20;
        REAL tol2=0.01;
        int numiter2=10;
        REAL l=0.5;
        REAL lambda0=1.;

        REAL FS  = IterativeProcessArcLength ( tol,numiter,tol2,numiter2,l,lambda0,conv );
        //REAL FS  = ShearRedNoIntegrationPointsArcLength ( numiter );

        return FS;
}

REAL SlopeAnalysis::SolveDeterministic ( bool IsSRM )
{
        REAL FSOLD,FS;
        int neqold;
        int neq=fCompMesh->NEquations();
        cout << "NUMBER OF EQUATIONS  = " << neq << endl;
        InitializeMemory();
        bool conv;

        if ( IsSRM==true ) {
                //FS = ShearRed ( 20,0.5,0.01 );
                FS =ShearRedNoIntegrationPoints( 20,0.5,0.01 );
        } else {
                //FS = GravityIncrease() ;

                FS  =ArcLength(conv);
                if(conv==false)
                {
                 //  FS = GravityIncrease() ;
                }
        }

        cout << "Refining.."<<endl;
        std::set<long> elindices,elindices2;
        for ( int iref=1; iref<=1; iref++ ) {
                cout << "computing deformation..."  << endl;
                ComputeElementDeformation();
                cout << "p refining..."  << endl;
                PRefineElementsAbove ( 0.001, fCompMesh->GetDefaultOrder()+iref,elindices2 );
                cout << "h refining..."  << endl;
                DivideElementsAbove ( 0.001,elindices );
                neqold=neq;
                neq=fCompMesh->NEquations();
                cout << "initializing memory..."  << endl;
                InitializeMemory();
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
                //   FS = GravityIncrease() ;
                }
                }
        }

        InitializeMemory();
        string meshref = "refinidemesh-grid";
        meshref+=".vtk";
        std::ofstream files ( meshref );
        TPZVTKGeoMesh::PrintGMeshVTK ( fCompMesh->Reference(),files,true );

        return FS;
}

REAL SlopeAnalysis::SolveSingleField(TPZVec<TPZFMatrix<REAL>> sample )
{
       REAL FSOLD;
        REAL tolfs=0.01;
        REAL refsqrt=0.01;
        bool conv;
        int neq=fCompMesh->NEquations();
        cout << "Starting to solve field  Mesh with "<< neq << " equations "<< endl;
        TransferFieldsSolutionFrom ( sample );

        //REAL FS = ShearRed ( 20,0.5,tolfs );
        REAL FS =ShearRedNoIntegrationPoints( 20,0.5,0.01 );
        //REAL FS  =ArcLength(conv);

        std::set<long> elindices,elindices2;
        for ( int iref=1; iref<=1; iref++ ) {
                cout << "refining level "<< iref <<endl;
                cout << "computing deformation..."  << endl;
                ComputeElementDeformation();
                cout << "p refining..."  << endl;
                PRefineElementsAbove ( refsqrt, fCompMesh->GetDefaultOrder()+iref,elindices2 );
                cout << "h refining..."  << endl;
                DivideElementsAbove ( refsqrt,elindices );
                neq=fCompMesh->NEquations();
                cout << "transfering solution..."  << endl;
                TransferFieldsSolutionFrom ( sample );
                cout <<  " Mesh with "<< neq << " equations "<< " fabs(FS-FSOLD)  "  << fabs ( FS-FSOLD )  << endl;
                FSOLD=FS;
               // FS = ShearRed ( 20,FSOLD,tolfs );
                FS =ShearRedNoIntegrationPoints( 20,FSOLD,0.01 );
                //REAL FS  =ArcLength(conv);

        }

        TransferFieldsSolutionFrom ( sample );
//         auto var=to_string ( ifield );
//         string meshref = "post/refinidemesh-grid";
//         meshref+=var;
//         meshref+=".vtk";
//         std::ofstream files ( meshref );
//         TPZVTKGeoMesh::PrintGMeshVTK ( fCompMesh->Reference(),files,true );

        return FS;

}
REAL SlopeAnalysis::SolveSingleField ( int ifield )
{
        REAL FSOLD;
        REAL tolfs=0.01;
        REAL refsqrt=0.01;
        bool conv;
        int neq=fCompMesh->NEquations();
        cout << "Starting to solve field   " << ifield << " Mesh with "<< neq << " equations "<< endl;
        TransferFieldsSolutionFrom ( ifield );

       // REAL FS = ShearRed ( 20,0.5,tolfs );
        REAL FS =ShearRedNoIntegrationPoints( 20,0.5,0.01 );
        //REAL FS  =ArcLength(conv);

        std::set<long> elindices,elindices2;
        for ( int iref=1; iref<=1; iref++ ) {
                cout << "refining level "<< iref <<endl;
                cout << "computing deformation..."  << endl;
                ComputeElementDeformation();
                cout << "p refining..."  << endl;
                PRefineElementsAbove ( refsqrt, fCompMesh->GetDefaultOrder()+iref,elindices2 );
                cout << "h refining..."  << endl;
                DivideElementsAbove ( refsqrt,elindices );
                neq=fCompMesh->NEquations();
                cout << "transfering solution..."  << endl;
                TransferFieldsSolutionFrom ( ifield );
                cout <<  " Mesh with "<< neq << " equations "<< " fabs(FS-FSOLD)  "  << fabs ( FS-FSOLD )  << endl;
                FSOLD=FS;
                //FS = ShearRed ( 20,FSOLD,tolfs );
                FS =ShearRedNoIntegrationPoints( 20,FSOLD,0.01 );
                //REAL FS  =ArcLength(conv);

        }

        TransferFieldsSolutionFrom ( ifield );
//         auto var=to_string ( ifield );
//         string meshref = "post/refinidemesh-grid";
//         meshref+=var;
//         meshref+=".vtk";
//         std::ofstream files ( meshref );
//         TPZVTKGeoMesh::PrintGMeshVTK ( fCompMesh->Reference(),files,true );

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

void SlopeAnalysis::SetFieldsData ( TPZCompMesh *CompMeshField, TPZFMatrix<REAL> SolutionValVec, TPZVec<REAL> meanvec, TPZVec<REAL> covvec, int samples )
{


        fSolutionValVec = SolutionValVec;
        fMeanvec = meanvec;
        fCovvec = covvec;
        fNSamples = samples;
        if ( fCompMeshField == CompMeshField ) {
                return;
        } else {
                fCompMeshField = CompMeshField;
        }
}
void SlopeAnalysis::SetFields ( TPZVec<TPZFMatrix<REAL>> fields )
{
        fFields = fields;
}

void SlopeAnalysis::SetFieldsSamples ( TPZVec<TPZFMatrix<REAL>> fields )
{
        fFieldSamples = fields;
}

TPZVec<TPZFMatrix<REAL>> SlopeAnalysis::GetFieldsSamples()
{
        return fFieldSamples;
}

TPZVec<TPZFMatrix<REAL>> SlopeAnalysis::GetFields()
{
        return fFields;
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
                conv  =anal.IterativeProcess2 ( cout,tol2, NumIter,  linesearch,  checkconv,iters );

//                 int numit1=20,numit2=10;
//                 REAL tolfs=1.e-2,tolrhs=1.e-3,l=0.5;
//                 FS = IterativeProcessArcLength ( tolfs,numit1,tolrhs,numit2,l,FS,conv );


                auto t2 = chrono::high_resolution_clock::now();
                auto ms_int = chrono::duration_cast<chrono::milliseconds> ( t2 - t1 );
                norm = Norm ( anal.Rhs() );
                //cout << "| step = " << counterout << " FS = "<< FS <<" tempo  iterproc = "<<ms_int.count() << " ms " << " conv?" << conv << " iters = " <<iters<< endl;


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
        //LoadingRamp ( 1.77 );
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
                FindRoot(conv);
                //if(conv==false)
                //{
                 //       conv  =anal.IterativeProcess2 ( cout,tol2, NumIter,  linesearch,  checkconv,iters );
                //}

//                  int numit1=20,numit2=10;
//                  REAL tolfs=1.e-2,tolrhs=1.e-3,l=0.5;
//                  FS = IterativeProcessArcLength ( tolfs,numit1,tolrhs,numit2,l,FS,conv );


                auto t2 = chrono::high_resolution_clock::now();
                auto ms_int = chrono::duration_cast<chrono::milliseconds> ( t2 - t1 );
                norm = Norm ( anal.Rhs() );
                if(counterout%5==0)
                {
                        cout << "| step = " << counterout << " FS = "<< FS <<" tempo  iterproc = "<<ms_int.count() << " ms " << " conv?" << conv << " iters = " <<iters<< endl;
                }


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
                bool conv  =anal.IterativeProcess2 ( cout,tol2, NumIter,  linesearch,  checkconv,iters );
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


//Tranfere a solucao nodal da malha cmesh para os pontos de integracao da malha fCompMesh. Este metodo e usado para transferir a solucao dos
// //campos estocasticos. O metodo findelement e caro, e custa muito ao monte carlo.
void SlopeAnalysis::TransferFieldsSolutionFrom ( TPZVec<TPZFMatrix<REAL>> sample )
{

        TPZMatWithMem<TPZElastoPlasticMem> *pMatWithMem2 = dynamic_cast<TPZMatWithMem<TPZElastoPlasticMem> *> ( fCompMesh->MaterialVec() [1] );
        if ( pMatWithMem2 ) {
                pMatWithMem2->SetUpdateMem ( true );
        } else {
                DebugStop();
        }

        TPZVec<TPZFMatrix<REAL>> fields(2);


        int nfields = sample.size();

        fields.Resize(nfields);

        TPZVec<TPZFMatrix<REAL>> samples ( nfields );

        for ( int ifield=0; ifield<nfields; ifield++ ) {
                fields[ifield] = GenerateRandomField ( fMeanvec[ifield],fCovvec[ifield],fSolutionValVec,sample[ifield] );
        }

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
                                REAL datarealvalue=datafield.sol[0][0];

                                mem.m_elastoplastic_state.fmatpropinit[imesh] = datarealvalue;
                                mem.m_elastoplastic_state.fmatprop[imesh] = datarealvalue;

                        }
                        pMatWithMem2->SetUpdateMem ( false );

                }

        }
//fCompMeshField->Solution().Redim(0,0);
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


void SlopeAnalysis::InitializeMemory ( )
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

                        mem.m_elastoplastic_state.fmatpropinit[0] = fCohesion;
                        mem.m_elastoplastic_state.fmatpropinit[1] = fAtrito;
                        //asscociativo
                        mem.m_elastoplastic_state.fmatpropinit[2] = fAtrito;

                        mem.m_elastoplastic_state.fmatprop[0] = fCohesion;
                        mem.m_elastoplastic_state.fmatprop[1] = fAtrito;
                        //asscociativo
                        mem.m_elastoplastic_state.fmatprop[2] = fAtrito;

                }

        }

        pMatWithMem2->SetUpdateMem ( false );

}

void SlopeAnalysis::IntegrateFieldOverARegion ( int imc )
{
        string saida = "postx/regionmean";
        auto var=to_string ( imc );
        saida+=var;
        saida+=".dat";
        ofstream out ( saida );

        //campos estocasticos
        TPZVec<TPZFMatrix<REAL>> fields=GetFields();

        int nfields = fields.size();

        TPZVec<REAL> integrationmeanvalues ( nfields );

        TPZVec<REAL> integrationvariance ( nfields );

        for ( int imesh=0; imesh<nfields; imesh++ ) {

                TPZFMatrix<REAL> solu ( fields[imesh].Rows(),1 );
                //REAL mean=0;
                int ndegreesoffredoom=fields[imesh].Rows();
                std::vector<double> v ( ndegreesoffredoom );
                for ( int inodalsol=0; inodalsol< ndegreesoffredoom; inodalsol++ ) {
                        solu ( inodalsol,0 ) =fields[imesh] ( inodalsol,imc );
                        v[inodalsol]=solu ( inodalsol,0 );
                        //mean+=solu(inodalsol,0);

                }

                if ( true ) {                              //calcula e imprime media e cov
                        double sum = std::accumulate ( v.begin(), v.end(), 0.0 );
                        double mean = sum / v.size();

                        double sq_sum = std::inner_product ( v.begin(), v.end(), v.begin(), 0.0 );
                        double stdev = std::sqrt ( sq_sum / v.size() - mean * mean );

                        cout << "MEDIA = "<< mean << endl;
                        cout << "COV = "<< stdev/mean << endl;
                }

                fCompMeshField->LoadSolution ( solu );

                //num elementos malha elastoplastica
                int nels =  fCompMeshField->NElements();
                REAL val=0.;
                int cout=0;
                std::vector<double> values;
                for ( int iel=0; iel<nels; iel++ ) {

                        TPZCompEl *cel = fCompMeshField->ElementVec() [iel];
                        TPZGeoEl * gel=cel->Reference();
                        TPZInterpolationSpace *intel = dynamic_cast<TPZInterpolationSpace *> ( cel );
                        if ( !cel || !intel ) {
                                continue;
                        }

                        TPZManVector<REAL,3> point ( 2,0. );

                        TPZMaterialDataT<REAL> data;
                        intel->InitMaterialData ( data );
                        data.fNeedsSol = true;
                        intel->ComputeRequiredData ( data,point );
                        REAL x=data.x[0];
                        REAL y=data.x[1] ;
                        REAL xc=42.;
                        REAL yc=44.;
                        REAL r1=13.5;
                        REAL r2=14.5;
                        REAL residual = ( x-xc ) * ( x-xc )+ ( y-yc ) * ( y-yc )-r2*r2;
                        REAL residual2 = ( x-xc ) * ( x-xc )+ ( y-yc ) * ( y-yc )-r1*r1;
                        TPZVec<REAL> qsi ( 2,0. ), sol;
                        REAL elsol=0.;
                        //if((38.<x<42.) && (28.<y<30.))
                        REAL solui=0.;
                        REAL area=0.;
                        if ( residual2>0&&residual<0 ) {

                                TPZIntPoints &rule = intel->GetIntegrationRule();
                                int np = rule.NPoints();
                                for ( int ip = 0; ip<np; ip++ ) {
                                        TPZManVector<REAL> point ( 2,0. );
                                        REAL weight;
                                        rule.Point ( ip, point, weight );
                                        intel->ComputeSolution ( point,data,false );
                                        weight*=fabs ( data.detjac );
                                        solui+=weight*data.sol[0][0] ;
                                        TPZFMatrix<REAL> jac,jacinv;
                                        TPZFMatrix<REAL> axes;
                                        REAL detjac;
                                        gel->Jacobian ( point, jac, axes, detjac, jacinv );
                                        area += weight;
                                }
                                values.push_back ( solui/area );
                        }
                }
                double sum = std::accumulate ( values.begin(), values.end(), 0.0 );
                double mean = sum / values.size();
                out <<  mean <<endl;
        }

}
void SlopeAnalysis::IntegrateFieldOverARegionB ( REAL refineaboveval,int imc )
{
        string saida = "postx2/regionmean";
        auto var=to_string ( imc );
        saida+=var;
        saida+=".dat";
        ofstream out ( saida );
        TPZManVector<REAL,3> findel ( 3,0. ),qsi ( 2,0. );


        TransferFieldsSolutionFrom ( imc );

        long nelem = fCompMesh->NElements();
        std::vector<double> valuescoes,valuesphi,veccriterio;
        for ( long iel=0; iel<nelem; iel++ ) {
                TPZCompEl *cel = fCompMesh->ElementVec() [iel];
                if ( !cel ) {
                        continue;
                }
                TPZGeoEl * gel = cel->Reference();

                TPZInterpolationSpace *intel = dynamic_cast<TPZInterpolationSpace *> ( cel );
                if ( !intel ) {
                        DebugStop();
                }
                TPZMatWithMem<TPZElastoPlasticMem> *pMatWithMem2 = dynamic_cast<TPZMatWithMem<TPZElastoPlasticMem> *> ( cel->Material() );
                if ( !pMatWithMem2 ) {
                        continue;
                }
                const TPZFMatrix<STATE> &elsol = fCompMesh->ElementSolution();

                // cout << "elsol.Get ( iel,0 )  = " <<  elsol.Get ( iel,0 )  << endl;
                //if ( elsol.Get ( iel,0 ) <=refineaboveval ) {
                //       continue;
                // }

                const TPZIntPoints &intpoints = intel->GetIntegrationRule();
                int nint = intpoints.NPoints();
                TPZManVector<REAL,3> point ( 2,0. );
                // TPZVec<REAL> point ( 3,0. );


                TPZMaterialDataT<REAL> data;
                intel->InitMaterialData ( data );
                data.fNeedsSol = true;
                intel->ComputeRequiredData ( data,point );
                //if(26<data.x[0]&&27<data.x[1]&& data.x[0]<45)
                //{

                // cout << "x "<< data.x[0]<< endl;
                // cout << "y "<< data.x[1]<< endl;

//                         ECoesion=23,
//         EAtrito=24,
//EStrainPlasticJ2    = 19,
                TPZVec<REAL> solc  = intel->IntegrateSolution ( 23 );
                TPZVec<REAL> solp  = intel->IntegrateSolution ( 24 );
                //TPZVec<REAL> solsqrtj2  = intel->IntegrateSolution(19);


                TPZIntPoints &rule = intel->GetIntegrationRule();
                int np = rule.NPoints();
                REAL area=0;
                REAL coes=0.;
                REAL phi=0.;
                for ( int ip = 0; ip<np; ip++ ) {
                        TPZManVector<REAL> point ( 2,0. );
                        REAL weight;
                        rule.Point ( ip, point, weight );
                        intel->ComputeSolution ( point,data,false );
                        weight*=fabs ( data.detjac );
                        //solui+=weight*data.sol[0][0] ;
                        TPZFMatrix<REAL> jac,jacinv;
                        TPZFMatrix<REAL> axes;
                        REAL detjac;
                        gel->Jacobian ( point, jac, axes, detjac, jacinv );
                        area += weight;
                        /*
                                                int indexplastic =data.intGlobPtIndex;
                                                TPZElastoPlasticMem &mem = pMatWithMem2->MemItem ( indexplastic );


                                                coes+=weight*mem.m_elastoplastic_state.fmatprop[0];
                                                phi+=weight*mem.m_elastoplastic_state.fmatprop[1];*/

                }
                REAL x=data.x[0];
                REAL y=data.x[1] ;
                REAL xc=43.;
                REAL yc=46.;
                REAL r1=15;
                REAL r2=20;
                REAL residual = ( x-xc ) * ( x-xc )+ ( y-yc ) * ( y-yc )-r2*r2;
                REAL residual2 = ( x-xc ) * ( x-xc )+ ( y-yc ) * ( y-yc )-r1*r1;
                TPZVec<REAL> qsi ( 2,0. ), sol;
                if ( residual2>0&&residual<0 ) {
                        REAL c=solc[0]/area;
                        REAL phi= solp[0]/area ;
                        valuescoes.push_back ( c );
                        valuesphi.push_back ( phi ) ;
                        veccriterio.push_back ( 2*c*cos ( phi ) ) ;
                }
//}
        }
        double sum = std::accumulate ( valuescoes.begin(), valuescoes.end(), 0.0 );
        double meancoes = sum / valuescoes.size();
        sum = std::accumulate ( valuesphi.begin(), valuesphi.end(), 0.0 );
        double meanphi = sum / valuesphi.size();
        sum = std::accumulate ( veccriterio.begin(), veccriterio.end(), 0.0 );
        double meancri = sum / veccriterio.size();
        out <<  meancoes <<endl;
        out <<  meanphi <<endl;
        out <<  meancri <<endl;



}


void SlopeAnalysis::ManageFieldCretion()
{
        int nfields = fMeanvec.size();

        fFields.resize ( nfields );

        TPZVec<TPZFMatrix<REAL>> samples ( nfields );

        for ( int ifield=0; ifield<nfields; ifield++ ) {
                TPZFMatrix<REAL>sample = CreateNormalStandardSamples();
                samples[ifield]=sample;
                fFields[ifield] = GenerateRandomField ( fMeanvec[ifield],fCovvec[ifield],fSolutionValVec,sample );

        }
        SetFieldsSamples ( samples );
}

void SlopeAnalysis::ManageFieldCretionSubSet()
{
        int nfields = fMeanvec.size();

        fFields.resize ( nfields );

        TPZVec<TPZFMatrix<REAL>> samples ( nfields );

        for ( int ifield=0; ifield<nfields; ifield++ ) {
                TPZFMatrix<REAL>sample = CreateNormalStandardSamples();
                samples[ifield]=sample;
                fFields[ifield] = GenerateRandomField ( fMeanvec[ifield],fCovvec[ifield],fSolutionValVec,sample );

        }
        SetFieldsSamples ( samples );
}


// void SlopeAnalysis::ManageFieldCretion()
// {
//         int nfields = fMeanvec.size();
//         fFields.resize ( nfields );
//         if ( !nfields ) DebugStop();
//         if ( fFieldSamples.size() ==0 ) DebugStop();
//         for ( int ifield=0; ifield<nfields; ifield++ ) {
//                 fFields[ifield] = GenerateRandomField ( fMeanvec[ifield],fCovvec[ifield],fSolutionValVec,fFieldSamples[ifield] );
//
//         }
//
// }



void SlopeAnalysis:: ManageFieldCretion ( std::vector<int>  fieldindexes )
{
        int nfields = fMeanvec.size();
        fFields.resize ( nfields );
        //fPesos.resize ( nfields );
        //fFieldSamples.Resize ( nfields );
        if ( !nfields ) DebugStop();
        int ndofs = fSolutionValVec.Rows();
        int chopedcollums=fieldindexes.size();
        TPZFMatrix<REAL> SolutionValVecSelected ( ndofs,chopedcollums );
        SolutionValVecSelected.Zero();
        for ( int iM=0; iM<chopedcollums; iM++ ) {
                for ( int idof=0; idof<ndofs; idof++ ) {
                        SolutionValVecSelected ( idof,iM ) =fSolutionValVec ( idof,fieldindexes[iM] );
                }
        }

        for ( int ifield=0; ifield<nfields; ifield++ ) {
                fFields[ifield] = GenerateRandomField ( fMeanvec[ifield],fCovvec[ifield],SolutionValVecSelected,fFieldSamples[ifield] );
        }
//cout << "sdasssss"<<endl;
}

void SlopeAnalysis:: ManageFieldCretion ( std::vector<std::vector<int>>  fieldindexes )
{
        int nfields = fMeanvec.size();
        fFields.resize ( nfields );
        //fPesos.resize ( nfields );
        //fFieldSamples.Resize ( nfields );
        if ( !nfields ) DebugStop();
        int ndofs = fSolutionValVec.Rows();



        for ( int ifield=0; ifield<nfields; ifield++ ) {
                int chopedcollums=fieldindexes[ifield].size();
                TPZFMatrix<REAL> SolutionValVecSelected ( ndofs,chopedcollums );
                SolutionValVecSelected.Zero();
                for ( int iM=0; iM<chopedcollums; iM++ ) {
                        for ( int idof=0; idof<ndofs; idof++ ) {
                                SolutionValVecSelected ( idof,iM ) =fSolutionValVec ( idof,fieldindexes[ifield][iM] );
                        }
                }


                fFields[ifield] = GenerateRandomField ( fMeanvec[ifield],fCovvec[ifield],SolutionValVecSelected,fFieldSamples[ifield] );
        }
//cout << "sdasssss"<<endl;
}


TPZFMatrix<REAL>  SlopeAnalysis::GenerateRandomField ( REAL mean, REAL cov,TPZFMatrix<REAL> valvec, TPZFMatrix<REAL> stdnormalsamples )
{
        //std::cout <<" mean = "<< mean  << " cov = " << cov<<std::endl;
        if ( valvec.Rows() <=1 ) {
                std::cout <<" no fSolutionValVec"<< std::endl;
                DebugStop();
        }

        //cout << "kkkkkksdas"<<endl;
        TPZFMatrix<REAL> hhat,hhat2;
        valvec.Multiply ( stdnormalsamples, hhat );
        REAL xi = sqrt ( log ( 1. + cov*cov ) );
        REAL lambda = log ( mean ) - 0.5*xi * xi;
        int M=valvec.Cols();

        //cout << "sdas"<<endl;
        for ( int i = 0; i < hhat.Rows(); i++ ) {           //ndof
                for ( int j = 0; j < hhat.Cols(); j++ ) {   //samples
                        hhat ( i,j ) = exp ( lambda + xi * hhat ( i,j ) );
                        REAL val=0.;
                        for ( int iM=0; iM<M; iM++ ) {
                                //   val+=stdnormalsamples(iM,j)*valvec(i,iM);
                        }
                        //hhat2 ( i,j ) = exp ( lambda + xi * val );
                }
        }

        return hhat;
}

TPZFMatrix<REAL> SlopeAnalysis::CreateNormalStandardSamples( )
{
        int M = fSolutionValVec.Cols();
        std::normal_distribution<double> distribution ( 0., 1. );

        TPZFMatrix<REAL>  THETA ( M, fNSamples, 0. );
        for ( int isample = 0; isample < fNSamples; isample++ ) {
                for ( int irdvar = 0; irdvar < M; irdvar++ ) {
                        std::random_device rd{};
                        std::mt19937 generator{ rd() };
                        REAL xic = distribution ( generator );
                        REAL xiphi = distribution ( generator );
                        THETA ( irdvar,isample ) = xic;
                }
        }
        return THETA;
}

std::vector<std::pair<int, double>>SlopeAnalysis::CrudeMonteCarlo(int a,int b)
{
        std::vector<std::pair<int, double>> fsvec;
        for(int imc=a;imc<b;imc++)
        {
                REAL fs = SolveSingleField ( imc );
                fsvec.emplace_back(imc, fs);
        }
        return fsvec;
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
        case 1: {
                //cout << "Solver called with TPZPardisoSolver\n";
                TPZSSpStructMatrix<STATE> SSpStructMatrix ( fCompMesh );
                SSpStructMatrix.SetNumThreads ( fNumThreads );
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
        REAL mc_cohesion    = coes;                         //kPa
        REAL mc_phi         = atrito;
        REAL mc_psi         = mc_phi;

        //elasticresponse.Print(std::cout);

        plasticmorh mohrcoulombplasticstep;

        mohrcoulombplasticstep.fYC.SetUp ( mc_phi, mc_psi, mc_cohesion, elasticresponse );

        mohrcoulombplasticstep.fER = elasticresponse;

        mohrcoulombplasticstep.SetStrengthReductionFactor(1.);

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

        material->SetWhichLoadVector ( 0 );                 //option to compute the total internal force vecor fi=(Bt sigma+ N (b+gradu))

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


void SlopeAnalysis::IntegrateFieldOverARegion ( REAL refineaboveval,int imc )
{
        string saida = "postx2/regionmean";
        auto var=to_string ( imc );
        saida+=var;
        saida+=".dat";
        ofstream out ( saida );
        TPZManVector<REAL,3> findel ( 3,0. ),qsi ( 2,0. );


        long nelem = fCompMesh->NElements();
        std::vector<double> valuescoes,valuesphi,vecsqrtj2;
        for ( long iel=0; iel<nelem; iel++ ) {
                TPZCompEl *cel = fCompMesh->ElementVec() [iel];
                if ( !cel ) {
                        continue;
                }
                TPZGeoEl * gel = cel->Reference();

                TPZInterpolationSpace *intel = dynamic_cast<TPZInterpolationSpace *> ( cel );
                if ( !intel ) {
                        DebugStop();
                }
                TPZMatWithMem<TPZElastoPlasticMem> *pMatWithMem2 = dynamic_cast<TPZMatWithMem<TPZElastoPlasticMem> *> ( cel->Material() );
                if ( !pMatWithMem2 ) {
                        continue;
                }
                const TPZFMatrix<STATE> &elsol = fCompMesh->ElementSolution();

                // cout << "elsol.Get ( iel,0 )  = " <<  elsol.Get ( iel,0 )  << endl;
                //if ( elsol.Get ( iel,0 ) <=refineaboveval ) {
                //       continue;
                // }

                const TPZIntPoints &intpoints = intel->GetIntegrationRule();
                int nint = intpoints.NPoints();
                TPZManVector<REAL,3> point ( 2,0. );
                // TPZVec<REAL> point ( 3,0. );


                TPZMaterialDataT<REAL> data;
                intel->InitMaterialData ( data );
                data.fNeedsSol = true;
                intel->ComputeRequiredData ( data,point );
                //if(26<data.x[0]&&27<data.x[1]&& data.x[0]<45)
                //{

                // cout << "x "<< data.x[0]<< endl;
                // cout << "y "<< data.x[1]<< endl;

//                         ECoesion=23,
//         EAtrito=24,
//EStrainPlasticJ2    = 19,
                TPZVec<REAL> solc  = intel->IntegrateSolution ( 23 );
                TPZVec<REAL> solp  = intel->IntegrateSolution ( 24 );
                TPZVec<REAL> solsqrtj2  = intel->IntegrateSolution ( 19 );


                TPZIntPoints &rule = intel->GetIntegrationRule();
                int np = rule.NPoints();
                REAL area=0;
                REAL coes=0.;
                REAL phi=0.;
                for ( int ip = 0; ip<np; ip++ ) {
                        TPZManVector<REAL> point ( 2,0. );
                        REAL weight;
                        rule.Point ( ip, point, weight );
                        intel->ComputeSolution ( point,data,false );
                        weight*=fabs ( data.detjac );
                        //solui+=weight*data.sol[0][0] ;
                        TPZFMatrix<REAL> jac,jacinv;
                        TPZFMatrix<REAL> axes;
                        REAL detjac;
                        gel->Jacobian ( point, jac, axes, detjac, jacinv );
                        area += weight;
                        /*
                                                int indexplastic =data.intGlobPtIndex;
                                                TPZElastoPlasticMem &mem = pMatWithMem2->MemItem ( indexplastic );


                                                coes+=weight*mem.m_elastoplastic_state.fmatprop[0];
                                                phi+=weight*mem.m_elastoplastic_state.fmatprop[1];*/

                }
                //cout << "mean sqrtj2 = "<< solsqrtj2[0]/area << endl;
                if ( solsqrtj2[0]/area >refineaboveval ) {
                        vecsqrtj2.push_back ( solsqrtj2[0]/area );
                        valuescoes.push_back ( solc[0]/area );
                        valuesphi.push_back ( solp[0]/area );
                }
//}
        }
        double sum = std::accumulate ( valuescoes.begin(), valuescoes.end(), 0.0 );
        double meancoes = sum / valuescoes.size();
        sum = std::accumulate ( valuesphi.begin(), valuesphi.end(), 0.0 );
        double meanphi = sum / valuesphi.size();
        sum = std::accumulate ( vecsqrtj2.begin(), vecsqrtj2.end(), 0.0 );
        double meanj2 = sum / valuesphi.size();
        out <<  meancoes <<endl;
        out <<  meanphi <<endl;
        out <<  meanj2 <<endl;


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
                        TPZStack<long> subsubels;
                        subintel->SetPreferredOrder ( porder );
                        subintel->Divide(subels[is],subsubels,0);
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

//                                 //mais um nivel de refinamento
//                                TPZStack<long> subsubels;
//                                subintel->Divide(subels[is],subsubels,0);
//                                cout << "subsubels.size() = "<< subsubels.size() <<endl;
//                                cout << "subels.size() = "<< subels.size() <<endl;
//                                  for ( int is2=0; is2<subsubels.size(); is2++ )
//                                  {
//                                          TPZCompEl *subsubcel = fCompMesh->ElementVec() [subsubels[is2]];
//                                          TPZInterpolationSpace *subsubintel = dynamic_cast<TPZInterpolationSpace *> ( subsubcel );
//                                          if ( !subintel )
//                                          {
//                                                  DebugStop();
//
//                                         }
//                                         subsubintel->SetPreferredOrder ( porder+1 );
//                                 }
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
                        REAL sqj2el = 0.;
                        REAL phivalplane=0.;
                        for ( int ind=0; ind<numind; ind++ ) {
                                int memoryindex = memindices[ind];
                                if ( memoryindex < 0 ) {
                                        continue;
                                }
                                TPZElastoPlasticMem &mem = pMatWithMem2->MemItem ( memindices[ind] );
                                TPZTensor<REAL> plastic =mem.m_elastoplastic_state.EpsP();
                                TPZTensor<REAL> total =mem.m_elastoplastic_state.EpsT();
                                TPZTensor<REAL> & Sigma = mem.m_sigma;
                                TPZTensor<REAL>::TPZDecomposed eigensystem;
                                Sigma.EigenSystem ( eigensystem );
                                REAL sig1=eigensystem.fEigenvalues[0];
                                REAL sig3=eigensystem.fEigenvalues[2];
                                REAL coes, atrito;
                                coes = mem.m_elastoplastic_state.fmatprop[0];
                                atrito = mem.m_elastoplastic_state.fmatprop[1];
                                REAL phi= ( sig1-sig3 )+ ( sig1+sig3 ) *sin ( atrito )-2*coes*cos ( atrito );
                                //sigma[0] - sigma[2] + (sigma[0] + sigma[2]) * sinphi - 2. * c*cosphi;
//                                 REAL phiyield =mem.m_phi;
                                if ( phi>=0 ) {
                                        //cout << "phiyield= "<< phi << endl;
                                        //sqj2el = phi;
                                }
                                //TPZVec<REAL> phi;
                                //pMatWithMem2->GetPlasticity().Phi(total,phi);
                                REAL J2 = plastic.J2();
                                REAL sqj2 = sqrt ( J2 );
                                REAL val=mem.m_elastoplastic_state.VolHardening();
                                sqj2el = max ( sqj2,sqj2el );
                                // phivalplane=phi[0];

                        }
                        fPlasticDeformSqJ2[el] = sqj2el;
                }
        }
        fCompMesh->SetElementSolution ( 0, fPlasticDeformSqJ2 );
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

        postprocdeter->PostProcess ( 1 );

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


TPZGeoMesh *  SlopeAnalysis::QuadGMesh ( int ref )
{
        const std::string name ( "Darcy Flow Slope" );

        TPZGeoMesh *gmesh  =  new TPZGeoMesh();

        gmesh->SetName ( name );

        gmesh->SetDimension ( 2 );

        TPZVec<REAL> coord ( 2 );

//         vector<vector<double>> co= {
//                 {0., 0.}, {75., 0.}, {75., 30.},{45., 30.},{35., 40.},{0.,40.},
//                 {35./3., 40.},{2 * 35/3., 40.},
//                 {30., 40.},{30., 30.}, {60.,30.},{2* 35./3.,2* 35/3.},
//                 {45., 2* 35/3.},{35./3., 35/3.}, {60., 35./3.}
//         };

        vector<vector<double>> co= {
                {0., 0.}, {70., 0.}, {70., 30.},{40., 30.},{30., 40.},{0.,40.},
                {10., 40.},{20., 40.},{25., 40.},
                {25., 30.}, {60.,30.},{20.,20.},
                {40., 20},{10, 10.}, {60., 10.}
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
        TopoLine[0] = 5;
        TopoLine[1] = 0;
        new TPZGeoElRefPattern< pzgeom::TPZGeoLinear> ( id, TopoLine, -5, *gmesh ); //left

        id++;
        TopoLine[0] = 3;
        TopoLine[1] = 4;
        new TPZGeoElRefPattern< pzgeom::TPZGeoLinear> ( id, TopoLine, -6, *gmesh ); //ramp

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


void SlopeAnalysis::Write ( TPZStream &buf, int withclassid ) const
{
        fSolutionValVec.Write ( buf,withclassid );
        fFieldSamples[0].Write ( buf,withclassid );
        fFieldSamples[1].Write ( buf,withclassid );
        fFields[0].Write ( buf,withclassid );
        fFields[1].Write ( buf,withclassid );

		int nc = fFieldSamplesSubSetFN.size();
        int sz=fFieldSamplesSubSetFN[0].size();

		buf.Write(&nc);
        buf.Write(&sz);

        cout << "nc write!!! = "<< nc <<endl;
        cout << "sz write!!! = "<< sz <<endl;

		for(int c=0; c<nc; c++)
        {
                for(int j=0;j<sz;j++)
                {
                        fFieldSamplesSubSetFN[c][j].Write(buf,withclassid);
                }
        }

}

void SlopeAnalysis::Read ( TPZStream &buf, void *context )
{

        fSolutionValVec.Read ( buf,context );
        fFieldSamples.resize ( 2 );
        fFieldSamples[0].Read ( buf,context );
        fFieldSamples[1].Read ( buf,context );
        fFields.resize ( 2 );
        fFields[0].Read ( buf,context );
        fFields[1].Read ( buf,context );

        int nc,sz;
		buf.Read(&nc);
        buf.Read(&sz);

        fFieldSamplesSubSetFN.resize(nc);
        cout << "nc = "<< nc <<endl;
        cout << "sz = "<< sz <<endl;
		for(int c=0; c<nc; c++)
        {
                fFieldSamplesSubSetFN[c].resize(sz);
                for(int j=0;j<sz;j++)
                {
                        fFieldSamplesSubSetFN[c][j].Read ( buf,context );
                }
        }

}

// 	template<class T>
// 	static void WriteObjects(TPZStream &buf, const TPZVec<T> &vec)
// 	{
// 		long c,nc = vec.NElements();
// 		buf.Write(&nc,1);
// 		for(c=0; c<nc; c++)
//             vec[c].Write(buf,0);
// 	}
//
// 	template<class T>
// 	static void WriteObjects(TPZStream &buf, const std::vector<T> &vec)
// 	{
// 		int c,nc = vec.size();
// 		buf.Write(&nc,1);
// 		for(c=0; c<nc; c++) vec[c].Write(buf,0);
// 	}
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
