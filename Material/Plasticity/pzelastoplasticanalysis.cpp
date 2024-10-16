//$Id: pzelastoplasticanalysis.cpp,v 1.27 2010-11-23 18:58:05 diogo Exp $
#include "pzelastoplasticanalysis.h"
#include "pzcmesh.h"
#include "pzvec.h"
#include "pzmanvector.h"
#include "checkconv.h"
#include "TPZMatElastoPlastic.h"
#include "tpzautopointer.h"
#include "pzcompelwithmem.h"
#include "TPZElastoPlasticMem.h"
#include "pzblockdiag.h"
#include "TPZSpStructMatrix.h"
#include "pzfstrmatrix.h"
#include "pzbdstrmatrix.h"
#include "pzstepsolver.h"
#include "TPZMaterial.h"
#include "TPZBndCondT.h"
#include "TPZMatElastoPlastic2D.h"

#include "pzbuildmultiphysicsmesh.h"

#include <map>
#include <set>
#include <stdio.h>
#include <fstream>

#include "TPZMatrixSolver.h"

#include "pzlog.h"

#ifdef PZ_LOG
static TPZLogger EPAnalysisLogger("pz.analysis.elastoplastic");
static TPZLogger loggertest("testing");
#endif

using namespace std;


TPZElastoPlasticAnalysis::TPZElastoPlasticAnalysis() : TPZLinearAnalysis(), fPrecond(NULL) {
	//Mesh()->Solution().Zero(); already performed in the nonlinearanalysis base class
	//fSolution.Zero();
}

TPZElastoPlasticAnalysis::TPZElastoPlasticAnalysis(TPZCompMesh *mesh,std::ostream &out) : TPZLinearAnalysis(mesh,true), fPrecond(NULL) {

	int numeq = fCompMesh->NEquations();
	fCumSol.Redim(numeq,1);
	fCumSol.Zero();
	fSolution.Redim(numeq,1);
	fSolution.Zero();

	LoadSolution();
}

TPZElastoPlasticAnalysis::~TPZElastoPlasticAnalysis()
{
	if(fPrecond)delete fPrecond;

#ifdef PZ_LOG
{
    if(EPAnalysisLogger.isDebugEnabled()){
        std::stringstream sout;
        sout << "<<< TPZElastoPlasticAnalysis::~TPZElastoPlasticAnalysis() *** Killing Object\n";
        LOGPZ_DEBUG(EPAnalysisLogger,sout.str().c_str());
    }
}
#endif
}

bool TPZElastoPlasticAnalysis::FindRoot(int & iters){

    REAL normrhs=10000,normdu=10000,normrhsn=10000,normdun=10000,normrhs0;
    TPZFMatrix<STATE> x(Solution()), dx(Solution());
    x.Zero();
    dx.Zero();
    REAL tol = 1.e-3;
    int n_it = 100;
    AssembleResidual();
    normrhs0=Norm(Rhs());

    for (int i = 1; i <= n_it; i++) {
        Assemble();
        Solve();
        if ( false)
        {
            TPZFMatrix<STATE> nextSol;
            REAL LineSearchTol = 0.001 * Norm ( fSolution );
            const int niter =30;
            this->LineSearch ( x, fSolution, nextSol, LineSearchTol, niter );
            x = nextSol;
            dx=fSolution;
        }else{
        dx = Solution();
        x += dx;
        }
        LoadSolution(x);

        normdun=normdu;
        normdu=Norm(dx);

        normrhsn=normrhs;
        normrhs = Norm(Rhs())/normrhs0;


        iters=i;
        if (normrhs<tol) {
//std::cout <<"iter = "<< i << " normrhs= " << normrhs<< " normrhsn= " << normrhsn<<" normdu= " << normdu<< " normdun= " << normdun<< std::endl;
            return true;
        }else if(i>4&&normrhsn<normrhs&&normdun<normdu){
           // std::cout << "Fail to converge. Divergent method." << std::endl;
            return false;
        }
    }

    //std::cout << " Not converged. Maximum number of iterations reached." << std::endl;
    return false;
}

REAL TPZElastoPlasticAnalysis::MyLineSearch(const TPZFMatrix<REAL> &Wn, const TPZFMatrix<REAL> &DeltaW, TPZFMatrix<REAL> &NextW, REAL RhsNormPrev, REAL &RhsNormResult, int niter, bool & converging)
{
    TPZFMatrix<REAL> Interval = DeltaW;
    REAL scalefactor = 1.;
    int iter = 0;
    do {
        Interval *= scalefactor;
        NextW = Wn;
        NextW += Interval;
        LoadSolution(NextW);
        AssembleResidual();
        RhsNormResult = Norm(fRhs);
        cout << " Normrhs  = " << RhsNormResult <<endl;
        scalefactor *= 0.5;
        iter++;
    } while (RhsNormResult > RhsNormPrev && iter < niter);
    if(fabs(RhsNormResult - RhsNormPrev)<1.e-6 )
    {
        converging=false;
    }
    else
    {
        converging=true;
    }
    scalefactor *= 2.;
	return scalefactor;
}

REAL TPZElastoPlasticAnalysis::LineSearch(const TPZFMatrix<REAL> &Wn, const TPZFMatrix<REAL> &DeltaW, TPZFMatrix<REAL> &NextW, REAL RhsNormPrev, REAL &RhsNormResult, int niter, bool & converging){

    TPZFMatrix<REAL> Interval = DeltaW;

#ifdef PZDEBUG
    {
        TPZLinearAnalysis::LoadSolution(Wn);
        AssembleResidual();
        STATE normprev = Norm(fRhs);
        if (fabs(normprev - RhsNormPrev) > 1.e-6) {
            std::stringstream sout;
            sout << "Norm of Wn " << Norm(Wn) << std::endl;
            sout << "Input previous norm " << RhsNormPrev << " Computed Norm " << normprev;
            LOGPZ_ERROR(EPAnalysisLogger, sout.str())
        }
    }
#endif
    REAL scalefactor = 1.;
    int iter = 0;
    do {
        Interval *= scalefactor;
        NextW = Wn;
        NextW += Interval;
        TPZLinearAnalysis::LoadSolution(NextW);
        AssembleResidual();
#ifdef PZDEBUGBIG
        {
            static int count = 0;
            {
                std::stringstream filename,varname;
                filename << "Sol." << count << ".txt";
                varname << "DelSol" << count << " = ";
                ofstream out(filename.str().c_str());
                Interval.Print(varname.str().c_str(),out,EMathematicaInput);
            }
            std::stringstream filename,varname;
            filename << "Rhs." << count << ".txt";
            varname << "Rhs" << count++ << " = ";
            ofstream out(filename.str().c_str());
            fRhs.Print(varname.str().c_str(),out,EMathematicaInput);
        }
#endif
        RhsNormResult = Norm(fRhs);
#ifndef PLASTICITY_CLEAN_OUT
        std::cout << "scale factor " << scalefactor << " residure norm " << RhsNormResult << std::endl;
#endif
        scalefactor *= 0.5;
        iter++;
    } while (RhsNormResult > RhsNormPrev && iter < niter);
    if(fabs(RhsNormResult - RhsNormPrev)<1.e-6 )
    {
        converging=false;
    }
    else
    {
        converging=true;
    }
    scalefactor *= 2.;
	return scalefactor;

}//void
/*
bool TPZElastoPlasticAnalysis::IterativeProcess2(std::ostream &out,REAL tol,int numiter, bool linesearch, bool checkconv) {

	int iter = 0;
	REAL error = 1.e10;
	int numeq = fCompMesh->NEquations();
	//Mesh()->Solution().Zero();
	//fSolution->Zero();



	TPZFMatrix<REAL> prevsol(fSolution);
	if(prevsol.Rows() != numeq) prevsol.Redim(numeq,1);

#ifdef PZ_LOG_keep
    {
        std::stringstream sout;
        fSolution.Print("Solution for checkconv",sout);
        LOGPZ_DEBUG(EPAnalysisLogger, sout.str())
    }
#endif

	if(checkconv){
		TPZVec<REAL> coefs(1,1.);
		TPZFMatrix<REAL> range(numeq,1,1.e-5);
		CheckConvergence(*this,fSolution,range,coefs);
	}

    Assemble();
    REAL RhsNormPrev = Norm(fRhs);
	bool linesearchconv=true;
	while(error > tol && iter < numiter) {

		//fSolution.Redim(0,0);
        REAL RhsNormResult = 0.;
        AssembleResidual();
		Solve();
		if (linesearch){
			TPZFMatrix<REAL> nextSol;
			const int niter = 30;
			this->LineSearch(prevsol, fSolution, nextSol, RhsNormPrev, RhsNormResult, niter,linesearchconv);
			fSolution = nextSol;
		}
		else{
			(TPZFMatrix<STATE> &)fSolution += prevsol;
            LoadSolution();
            AssembleResidual();
            RhsNormResult = Norm(fRhs);
		}

		prevsol -= fSolution;
		REAL normDeltaSol = Norm(prevsol);
		prevsol = fSolution;
		REAL norm = RhsNormResult;
        RhsNormPrev = RhsNormResult;
		//       out << "Iteracao n : " << (iter+1) << " : norma da solucao |Delta(Un)|: " << norm << endl;
        std::cout << "Iteracao n : " << (iter+1) << " : normas |Delta(Un)| e |Delta(rhs)| : " << normDeltaSol << " / " << RhsNormResult << endl;
        //        std::cout << "Iteracao n : " << (iter+1) << " : fRhs : " << fRhs << endl;

		if(norm < tol) {
            std::cout << "\nTolerancia atingida na iteracao : " << (iter+1) << endl;
            std::cout << "\n\nNorma da solucao |Delta(Un)|  : " << norm << endl << endl;
            return true;

		} else
			if( (norm - error) > 1.e-9 ) {
                std::cout << "\nDivergent Method \n";
                return false;

			}
		error = norm;
		iter++;
		out.flush();
	}

}*/

bool TPZElastoPlasticAnalysis::IterativeProcess ( std::ostream &out,REAL tol,int numiter, bool linesearch, bool checkconv,int &iters )
{

    int iter = 0;
    REAL errordisplace = 1.e10,errorrhs=1.e10;
    int numeq = fCompMesh->NEquations();

    cout << "number of equations = " << numeq <<endl;

    TPZFMatrix<STATE> prevsol ( fSolution );
    if ( prevsol.Rows() != numeq ) prevsol.Redim ( numeq,1 );

    if ( checkconv )
    {
        TPZVec<REAL> coefs ( 1,1. );
        TPZFMatrix<STATE> range ( numeq,1,1. );
        CheckConvergence ( *this,fSolution,range,coefs );
    }
    bool a=true,b=true,c=true;

    Assemble();

    REAL normrhs0 = Norm ( fRhs );
    cout << "normrhs0 = " << normrhs0 << endl;

    while ( a  && c )
    {
        if(iter%1==0)
        {
            //cout<< "Assembling in iter = "<<iter<< endl;
            Assemble();
        }
        Solve();
        if ( linesearch )
        {
            TPZFMatrix<STATE> nextSol;
            //REAL LineSearchTol = 1e-3 * Norm(fSolution);
            REAL LineSearchTol = 0.001 * Norm ( fSolution );
            const int niter =10;
            this->LineSearch ( prevsol, fSolution, nextSol, LineSearchTol, niter );
            fSolution = nextSol;
        }
        else
        {
            TPZFMatrix<STATE> sol = fSolution;
            sol += prevsol;
        }

        prevsol -= fSolution;
        REAL normu = Norm ( prevsol );

        prevsol = fSolution;
        this->LoadSolution ( fSolution );
        //this->AssembleResidual();

        REAL normf  =  Norm ( fRhs )/normrhs0;
        cout << "Iteracao n : " << ( iter ) << " : normas |Delta(Un)| e |Delta(rhs)/rhs0| : " << normu << " / " <<normf<< " | tol = "<<tol << endl;
        a = iter < numiter ;
        b =errordisplace > tol;
        c= errorrhs > tol;

        if ( ( iter >=numiter || ( iter>10&& normf >errorrhs &&normu>errordisplace) ) )
        {
            cout << "\nDivergent Method\n";
            return false;
        }
        errorrhs = normf;
        errordisplace=normu;
        iter++;
        out.flush();

    }
    iters=iter;
    cout << "Iteracao n : " << ( iter ) << "Norm ( prevsol ) = "<<Norm ( prevsol ) << "Norm ( fRhs ) = "<<Norm ( fRhs ) << endl;
    return true;
}
bool TPZElastoPlasticAnalysis::IterativeProcess2 ( std::ostream &out,REAL tol,int numiter, bool linesearch, bool checkconv,int &iters )
{

    int iter = 0;
    REAL errordisplace = 1.e10,errorrhs=1.e10;
    int numeq = fCompMesh->NEquations();

    //cout << "number of equations = " << numeq <<endl;

    TPZFMatrix<STATE> prevsol ( fSolution );
    if ( prevsol.Rows() != numeq ) prevsol.Redim ( numeq,1 );

    if ( checkconv )
    {
        TPZVec<REAL> coefs ( 1,1. );
        TPZFMatrix<STATE> range ( numeq,1,1. );
        CheckConvergence ( *this,fSolution,range,coefs );
    }
    bool a=true,b=true,c=true;

    Assemble();

    REAL normrhs0 = Norm ( fRhs );
    //cout << "normrhs0 = " << normrhs0 << endl;

    while ( a  && c &&b)
    {

        Assemble();
        Solve();

        if ( linesearch )
        {
            TPZFMatrix<STATE> nextSol;
            //REAL LineSearchTol = 1e-3 * Norm(fSolution);
            REAL LineSearchTol = 0.001 * Norm ( fSolution );
            const int niter =10;
            this->LineSearch ( prevsol, fSolution, nextSol, LineSearchTol, niter );
            fSolution = nextSol;
        }
        else
        {
            TPZFMatrix<STATE> sol = fSolution;
            sol += prevsol;
        }

        prevsol -= fSolution;
        REAL normu = Norm ( prevsol );

        prevsol = fSolution;
        this->LoadSolution ( fSolution );
        //this->AssembleResidual();

        REAL normf  =  Norm ( fRhs )/normrhs0;
        //cout << "Iteracao n : " << ( iter ) << " : normas |Delta(Un)| e |Delta(rhs)/rhs0| : " << normu << " / " <<normf<< " | tol = "<<tol << endl;
        a = iter < numiter ;
        b =errordisplace > tol;
        c= errorrhs > tol;

        iters=iter;
        if ( ( iter >=numiter || ( iter>2&& normf >errorrhs &&normu>errordisplace) ) )
        {
            //cout << "\nDivergent Method\n";
            return false;
        }
        errorrhs = normf;
        errordisplace=normu;
        iter++;
        out.flush();

    }

    //cout << "Iteracao n : " << ( iter ) << "Norm ( prevsol ) = "<<Norm ( prevsol ) << "Norm ( fRhs ) = "<<Norm ( fRhs ) << endl;
    return true;
}

void TPZElastoPlasticAnalysis::TransferSolution()
{


}

REAL TPZElastoPlasticAnalysis::LineSearch ( const TPZFMatrix<STATE> &Wn, TPZFMatrix<STATE> DeltaW, TPZFMatrix<STATE> &NextW, REAL tol, int niter )
{

    //cout << "Entering line search "<<endl;
    REAL error = 2.*tol+1.;
    REAL A = 0.1, B = 2., L = 0, M = 0.;
    TPZFMatrix<STATE> ak, bk, lambdak, muk, Interval;
    REAL NormResLambda = 0., NormResMu = 0.;
    //ak = Wn + 0.1 * DeltaW
    ak = DeltaW;
    ak *= A;
    ak += Wn;
    //bk = Wn + 2. DeltaW
    bk = DeltaW;
    bk *= B;
    bk += Wn;
    //Interval = (bk-ak)
    Interval = bk;
    Interval -= ak;
    int iter = 0;
    int KeptVal = -1; //0 means I have residual(labmda); 1 means I have residual(mu); -1 means I have nothing
    while ( error > tol && iter < niter )
    {
        iter++;
        //cout << "a  " << std::endl;
        if ( KeptVal != 0 )
        {
            L = 0.382* ( B-A )+A;
            //lambdak = ak + 0.382*(bk-ak)
            lambdak = Interval;
            lambdak *= 0.382;
            lambdak += ak;
            //computing residual
            LoadSolution ( lambdak );
            this->AssembleResidual();
            NormResLambda = Norm ( fRhs );
        }

        if ( KeptVal != 1 )
        {
            //muk = ak + 0.618*(bk-ak)
            M = 0.618* ( B-A )+A;
            muk = Interval;
            muk *= 0.618;
            muk += ak;
            LoadSolution ( muk );
            this->AssembleResidual();
            NormResMu = Norm ( fRhs );
        }

        if ( NormResLambda > NormResMu )
        {
            A = L;
            L = M;
            ak = lambdak;
            lambdak = muk;
            NormResLambda = NormResMu;
            KeptVal = 0;
        }
        else
        {
            B = M;
            M = L;
            bk = muk;
            muk = lambdak;
            NormResMu = NormResLambda;
            KeptVal = 1;
        }
        //error = Norm(bk-ak)
        Interval = bk;
        Interval -= ak;
        error = Norm ( Interval );

        //alpha shall be alpha <= 1
        if ( A > 1. && B > 1. ) break;

    }//while

    double ALPHA = 0.5* ( A + B );
    NextW = ak;
    NextW += bk;
    NextW *= 0.5;

    if ( ALPHA > 1. ) //alpha shall be alpha <= 1
    {
        NextW = Wn;
        NextW += DeltaW;
        return 1.;
    }

    return ALPHA;

}//void

void TPZElastoPlasticAnalysis::SetUpdateMem(int update)
{
	if(!fCompMesh)return;

	std::map<int, TPZMaterial *> & refMatVec = fCompMesh->MaterialVec();

    std::map<int, TPZMaterial * >::iterator mit;

	TPZMatWithMem<TPZElastoPlasticMem> * pMatWithMem; // defined in file pzelastoplastic.h
	TPZMatWithMem<TPZPoroElastoPlasticMem> * pMatWithMem2; // define in file pzporous.h

//    TPZMatElastoPlasticSest2D< TPZElasticCriteria >

    for(mit=refMatVec.begin(); mit!= refMatVec.end(); mit++)
    {
        pMatWithMem = dynamic_cast<TPZMatWithMem<TPZElastoPlasticMem> *>( mit->second );
		if(pMatWithMem != NULL)
        {
           pMatWithMem->SetUpdateMem(update);
        }
        pMatWithMem2 = dynamic_cast<TPZMatWithMem<TPZPoroElastoPlasticMem> *>( mit->second);
		if(pMatWithMem2 != NULL)
        {
            pMatWithMem2->SetUpdateMem(update);
        }
    }

}

#include "Elasticity/TPZElasticity2D.h"

REAL TPZElastoPlasticAnalysis::AcceptSolution(const int ResetOutputDisplacements)
{

    TPZMaterial *mat = fCompMesh->FindMaterial(1);
    if (!mat) {
        DebugStop();
    }
    auto *elasmat = dynamic_cast<TPZElasticity2D *>(mat);
    if(elasmat)
    {
        // the material is linear
        return 0.;
    }


	if(ResetOutputDisplacements)
	{
		fCumSol.Zero();
	}else{
		fCumSol += fSolution;
	}

	#ifdef PZ_LOG
	{
            if (EPAnalysisLogger.isDebugEnabled()){
               std::stringstream sout;
               sout << ">>> TTPZElastoPlasticAnalysis::AcceptSolution *** "
                    << " with Norm(fCumSol) = " << Norm(fCumSol);
               LOGPZ_DEBUG(EPAnalysisLogger,sout.str().c_str());
            }
	}
	#endif

	this->SetUpdateMem(true);

	fRhs.Zero();

    AssembleResidual();
	REAL norm = Norm(fRhs);

	this->SetUpdateMem(false);

	fSolution.Zero();

	LoadSolution();


	return norm;
}

/** @brief Load the solution into the computable grid, transferring it to the multi physics meshes */
void TPZElastoPlasticAnalysis::LoadSolution()
{
    TPZLinearAnalysis::LoadSolution();
    //a verificacao retorna verdadeiro ou falso para: return fMultiPhysics != NULL;
    //cout << this->IsMultiPhysicsConfiguration() << endl;
        if (this->IsMultiPhysicsConfiguration()) {
            cout << "nao é multifisica, porque entra aqui?" <<endl;
        //TPZBuildMultiphysicsMesh::TransferFromMultiPhysics(fMeshVec, fMultiPhysics);
            //fCompMesh->TransferMultiphysicsSolution();?
    }

}



void TPZElastoPlasticAnalysis::CheckConv(std::ostream &out, REAL range) {

#ifdef PZ_LOG
{
   std::stringstream sout;
   sout << ">>> TPZElastoPlasticAnalysis::CheckConv() ***"
        << "\nEntering method with parameters:"
	    << "\n range = " << range;
   LOGPZ_INFO(EPAnalysisLogger,sout.str().c_str());
}
#endif

   int numeq = fCompMesh->NEquations();

   TPZFMatrix<REAL> rangeMatrix(numeq, 1, range);

   TPZVec<REAL> coefs(1,1.);

   CheckConvergence(*this,fSolution,rangeMatrix,coefs);

}

void TPZElastoPlasticAnalysis::ComputeTangent(TPZFMatrix<REAL> &tangent, TPZVec<REAL> &coefs, int icase){

	int neq = fCompMesh->NEquations();
	tangent.Redim(neq,neq);
	TPZFMatrix<REAL> rhs(neq,1);
	TPZFStructMatrix<STATE> substitute(Mesh());
	TPZAutoPointer<TPZGuiInterface> guiInterface(0);
	substitute.Assemble(tangent,rhs,guiInterface);
//	TPZStructMatrix::Assemble(tangent, rhs, *Mesh());
}

int TPZElastoPlasticAnalysis::NumCases(){
	return 1;
}

void TPZElastoPlasticAnalysis::Residual(TPZFMatrix<REAL> &residual, int icase){
	int neq = fCompMesh->NEquations();
//	TPZFMatrix<REAL> tangent(neq,neq);
	residual.Redim(neq,1);
	TPZFStructMatrix<STATE> substitute(Mesh());
	TPZAutoPointer<TPZGuiInterface> guiInterface(0);
	substitute.Assemble(residual,guiInterface);
//	TPZStructMatrix::Assemble(/*tangent,*/ residual, *Mesh());
	residual *= -1;
}

void TPZElastoPlasticAnalysis::SetPrecond(TPZMatrixSolver<REAL> &precond){
  if(fPrecond) delete fPrecond;
    fPrecond = (TPZMatrixSolver<REAL> *) precond.Clone();
}

void TPZElastoPlasticAnalysis::UpdatePrecond()
{
   if(fPrecond)
   {
       TPZMatrix<REAL> * pMatrix = TPZLinearAnalysis::MatrixSolver<STATE>().Matrix().operator->();
		TPZMatrix<REAL> * pPrecondMat = fPrecond->Matrix().operator->();
		pPrecondMat->Zero();
		TPZBlockDiagonal<REAL> *pBlock = dynamic_cast<TPZBlockDiagonal<REAL> *>(pPrecondMat);
		pBlock->BuildFromMatrix(*pMatrix);
   }
}

void TPZElastoPlasticAnalysis::SetBiCGStab(int numiter, REAL tol)
{
#ifdef PZ_LOG
{
   std::stringstream sout;
   sout << ">>> TPZElastoPlasticAnalysis::SetBiCGStab() *** numiter = " << numiter << " and tol=" << tol;
   LOGPZ_INFO(EPAnalysisLogger,sout.str().c_str());
}
#endif

	TPZSpStructMatrix<STATE> StrMatrix(Mesh());
    this->SetStructuralMatrix(StrMatrix);
	TPZMatrix<REAL> * mat = StrMatrix.Create();

    TPZBlockDiagonalStructMatrix<STATE> strBlockDiag(Mesh());
    TPZStepSolver<REAL> Pre;
    TPZBlockDiagonal<REAL> * block = new TPZBlockDiagonal<REAL>();

#ifdef PZ_LOG
{
   std::stringstream sout;
   sout << "*** TPZElastoPlasticAnalysis::SetBiCGStab() *** Assembling Block Diagonal Preconditioning matrix\n";
   LOGPZ_INFO(EPAnalysisLogger,sout.str().c_str());
}
#endif

    strBlockDiag.AssembleBlockDiagonal(*block); // just to initialize structure
	Pre.SetMatrix(block);
    Pre.SetDirect(ELU);
    TPZStepSolver<REAL> Solver;
 	Solver.SetBiCGStab(numiter, Pre, tol, 0);
    Solver.SetMatrix(mat);
    this->SetSolver(Solver);
	this->SetPrecond(Pre);

#ifdef PZ_LOG
{
   std::stringstream sout;
   sout << "<<< TPZElastoPlasticAnalysis::SetBiCGStab() *** Exiting\n";
   LOGPZ_INFO(EPAnalysisLogger,sout.str().c_str());
}
#endif

}


void TPZElastoPlasticAnalysis::SetBiCGStab_Jacobi(int numiter, REAL tol)
{
#ifdef PZ_LOG
{
   std::stringstream sout;
   sout << ">>> TPZElastoPlasticAnalysis::SetBiCGStab_Jacobi() *** numiter = " << numiter << " and tol=" << tol;
   LOGPZ_INFO(EPAnalysisLogger,sout.str().c_str());
}
#endif

	TPZSpStructMatrix<STATE> StrMatrix(Mesh());
//	TPZFStructMatrix StrMatrix(Mesh());
    this->SetStructuralMatrix(StrMatrix);
	TPZMatrix<REAL> * mat = StrMatrix.Create();

    TPZBlockDiagonalStructMatrix<STATE> strBlockDiag(Mesh());
    TPZStepSolver<REAL> Pre;
    TPZBlockDiagonal<REAL> * block = new TPZBlockDiagonal<REAL>();

#ifdef PZ_LOG
{
   std::stringstream sout;
   sout << "*** TPZElastoPlasticAnalysis::SetBiCGStab_Jacobi() *** Assembling Block Diagonal Preconditioning matrix\n";
   LOGPZ_INFO(EPAnalysisLogger,sout.str().c_str());
}
#endif

    strBlockDiag.AssembleBlockDiagonal(*block); // just to initialize structure
	Pre.SetMatrix(block);
    //    Pre.SetDirect(ELU);
    //Pre.SetDirect(ELDLt);
	Pre.SetJacobi(numiter, tol, 0);
    TPZStepSolver<REAL> Solver;
 	Solver.SetBiCGStab(numiter, Pre, tol, 0);
    Solver.SetMatrix(mat);
    this->SetSolver(Solver);
	this->SetPrecond(Pre);

#ifdef PZ_LOG
{
   std::stringstream sout;
   sout << "<<< TPZElastoPlasticAnalysis::SetBiCGStab_Jacobi() *** Exiting\n";
   LOGPZ_INFO(EPAnalysisLogger,sout.str().c_str());
}
#endif
}

void TPZElastoPlasticAnalysis::SetLU()
{
#ifdef PZ_LOG
{
   std::stringstream sout;
   sout << ">>> TPZElastoPlasticAnalysis::SetLU() ***\n";
   LOGPZ_INFO(EPAnalysisLogger,sout.str().c_str());
}
#endif

    TPZFStructMatrix<STATE> StrMatrix(Mesh());
    this->SetStructuralMatrix(StrMatrix);

    TPZMatrix<REAL> * mat = StrMatrix.Create();

    TPZStepSolver<REAL> Solver;
    //Solver.SetDirect(ELU);// ECholesky -> simétrica e positiva definida
	Solver.SetDirect(ELU);
    Solver.SetMatrix(mat);

    this->SetSolver(Solver);
}

void TPZElastoPlasticAnalysis::TransferSolution(TPZPostProcAnalysis & ppanalysis)
{
	TPZFMatrix<REAL> bkpSolution = fSolution;


	fSolution = fCumSol;
//	 fSolution.Print();
	LoadSolution();//Carrega a solucao convergida no analysis
	//passa o cum sol para o post
	ppanalysis.TransferSolution();//Transfere solucao convergida para o pos processamento


	fSolution = bkpSolution;

	LoadSolution();
}

void TPZElastoPlasticAnalysis::ManageIterativeProcess(std::ostream &out,REAL tol,int numiter,
									int BCId, int nsteps, REAL PGRatio,
									TPZFMatrix<REAL> & val1Begin, TPZFMatrix<REAL> & val1End,
									TPZFMatrix<REAL> & val2Begin, TPZFMatrix<REAL> & val2End,
									TPZPostProcAnalysis * ppAnalysis, int res)
{

	if(!fCompMesh)return;

#ifdef PZ_LOG
{

   std::stringstream sout;
   sout << "<<< TPZElastoPlasticAnalysis::ManageIterativeProcess() ***";
   sout << "\nWith parameters:\n";
   sout << "\ntol = " << tol;
   sout << "\nnumiter = " << numiter;
   sout << "\nBCId = " << BCId;
   sout << "\nnsteps = " << nsteps;
   sout << "\nPGRatio = " << PGRatio;
   sout << "\nval1Begin = " << val1Begin;
   sout << "\nval1End = " << val1End;
   sout << "\nval2Begin = " << val2Begin;
   sout << "\nval2End = " << val2End;
   if(ppAnalysis)
	{
		sout << "\nppanalysis set";
	}else
	{
		sout << "\nppanalysis NOT set";
	}
   LOGPZ_INFO(EPAnalysisLogger,sout.str().c_str());
}
#endif

	// computing the initial value for the PG progression such that its sum equals one;
	REAL a0;
	if(fabs(PGRatio - 1.) < 1.e-3)
	{
	    a0 = 1. / REAL(nsteps);
	}else{
		a0 = (PGRatio - 1) / (pow(PGRatio,nsteps) - 1.);
	}
	TPZFNMatrix<36> val1(6,6,0.), deltaVal1(6,6,0.);
	TPZFNMatrix< 6> val2(6,1,0.), deltaVal2(6,1,0.);

	deltaVal1 = val1End;
	deltaVal1.ZAXPY(-1., val1Begin);
	deltaVal2 = val2End;
	deltaVal2.ZAXPY(-1., val2Begin);

	// ZAXPY operation: *this += alpha * p

	TPZMaterial * mat = fCompMesh->FindMaterial(BCId);
	auto * pBC = dynamic_cast<TPZBndCondT<STATE> *>(mat);
	if(!pBC)return;

    int i;
	for(i = 0; i < nsteps; i++)
	{
		REAL stepLen;
		if(fabs(PGRatio - 1.) < 1.e-3)
		{
			stepLen = REAL(i+1) / REAL(nsteps);
		}else{
		    stepLen = a0 * (pow(PGRatio,i+1) - 1) / (PGRatio - 1.);
		}

		val1 = val1Begin;
		val1.ZAXPY(stepLen, deltaVal1);
		val2 = val2Begin;
		val2.ZAXPY(stepLen, deltaVal2);
		TPZManVector<STATE,6> actualVal2(6,0);
        for(int i = 0; i < 6; i++) actualVal2[i]=val2(i,0);
		pBC->SetVal1(val1);
		pBC->SetVal2(actualVal2);

		#ifdef PZ_LOG
		{
		   std::stringstream sout;
		   sout << "*** TPZElastoPlasticAnalysis::ManageIterativeProcess() *** load step " << i;
		   sout << " stepLen = " << stepLen;
		   sout << "\nBC.val1() = " << val1;
		   sout << "\nBC.val2() = " << val2;
		   LOGPZ_INFO(EPAnalysisLogger,sout.str().c_str());
		}
		#endif

        bool linesearch = false;
        bool checkconv = false;
            bool convordiv;
            int iters;
		IterativeProcess(out, tol, numiter, linesearch, checkconv,iters);


		#ifdef PZ_LOG
		{
		   std::stringstream sout;
		   sout << "*** TPZElastoPlasticAnalysis::ManageIterativeProcess() *** load step " << i << " ended";
		   LOGPZ_INFO(EPAnalysisLogger,sout.str().c_str());
		}
		#endif

		AcceptSolution();

		if(ppAnalysis)
		{
			#ifdef PZ_LOG
			{
			   std::stringstream sout;
			   sout << "*** TPZElastoPlasticAnalysis::ManageIterativeProcess() *** PostProcessing ";
			   LOGPZ_INFO(EPAnalysisLogger,sout.str().c_str());
			}
			#endif
			TransferSolution(*ppAnalysis);
			ppAnalysis->PostProcess(res);
		}
	}

	#ifdef PZ_LOG
	{
	   std::stringstream sout;
	   sout << "<<< TPZElastoPlasticAnalysis::ManageIterativeProcess() *** Exiting";
	   LOGPZ_INFO(EPAnalysisLogger,sout.str().c_str());
	}
	#endif
}

// CompEl create Functions setup

#include "pzintel.h"
//#include "pzelctempplus.h"

#include "pzrefpoint.h"
#include "pzgeopoint.h"
#include "pzshapepoint.h"
#include "tpzpoint.h"

#include "pzshapelinear.h"
#include "TPZGeoLinear.h"
#include "TPZRefLinear.h"
#include "tpzline.h"

#include "pzshapetriang.h"
#include "pzreftriangle.h"
#include "pzgeotriangle.h"
#include "tpztriangle.h"

#include "pzrefquad.h"
#include "pzshapequad.h"
#include "pzgeoquad.h"
#include "tpzquadrilateral.h"

#include "pzshapeprism.h"
#include "pzrefprism.h"
#include "pzgeoprism.h"
#include "tpzprism.h"

#include "pzshapetetra.h"
#include "pzreftetrahedra.h"
#include "pzgeotetrahedra.h"
#include "tpztetrahedron.h"

#include "pzshapepiram.h"
#include "pzrefpyram.h"
#include "pzgeopyramid.h"
#include "tpzpyramid.h"

#include "TPZGeoCube.h"
#include "pzshapecube.h"
#include "TPZRefCube.h"
#include "tpzcube.h"

#include "pzelctemp.h"

#include "TPZCompElH1.h"
void TPZElastoPlasticAnalysis::SetAllCreateFunctionsWithMem(TPZCompMesh *cmesh)
{
/*	pzgeom::TPZGeoPoint::fp = TPZElastoPlasticAnalysis::CreatePointElWithMem;
	 pzgeom::TPZGeoQuad::fp = TPZElastoPlasticAnalysis::CreateQuadElWithMem;
	pzgeom::TPZGeoTriangle::fp = TPZElastoPlasticAnalysis::CreateTriangElWithMem;
	pzgeom::TPZGeoPrism::fp = TPZElastoPlasticAnalysis::CreatePrismElWithMem;
	pzgeom::TPZGeoTetrahedra::fp = TPZElastoPlasticAnalysis::CreateTetraElWithMem;
	pzgeom::TPZGeoPyramid::fp = TPZElastoPlasticAnalysis::CreatePyramElWithMem;
	pzgeom::TPZGeoCube::fp = TPZElastoPlasticAnalysis::CreateCubeElWithMem;
*/
 /*   TPZManVector<TCreateFunction,10> functions(8);
    functions[EPoint] = &TPZElastoPlasticAnalysis::CreatePointElWithMem;
	functions[EOned] = TPZElastoPlasticAnalysis::CreateLinearElWithMem;
	functions[EQuadrilateral] = TPZElastoPlasticAnalysis::CreateQuadElWithMem;
	functions[ETriangle] = TPZElastoPlasticAnalysis::CreateTriangElWithMem;
	functions[EPrisma] = TPZElastoPlasticAnalysis::CreatePrismElWithMem;
	functions[ETetraedro] = TPZElastoPlasticAnalysis::CreateTetraElWithMem;
	functions[EPiramide] = TPZElastoPlasticAnalysis::CreatePyramElWithMem;
	functions[ECube] = TPZElastoPlasticAnalysis::CreateCubeElWithMem;
    */
 TPZManVector<TCreateFunction,10> functions(8);
	TCreateFunction fp[8];
    cmesh->ApproxSpace().SetCreateFunctions(functions);

}

TPZCompEl * TPZElastoPlasticAnalysis::CreateCubeElWithMem(TPZGeoEl *gel, TPZCompMesh &mesh, int64_t &index)
{
	//TPZCompElWithMem<TPZCompElH1<pzshape::TPZShapeCube> >
	return new TPZCompElWithMem<TPZCompElH1<pzshape::TPZShapeCube> >(mesh,gel);
}

TPZCompEl * TPZElastoPlasticAnalysis::CreateLinearElWithMem(TPZGeoEl *gel, TPZCompMesh &mesh, int64_t &index)
{
	return new TPZCompElWithMem<TPZCompElH1<pzshape::TPZShapeLinear > >(mesh,gel);
}

TPZCompEl * TPZElastoPlasticAnalysis::CreatePointElWithMem(TPZGeoEl *gel, TPZCompMesh &mesh, int64_t &index)
{
	return new TPZCompElWithMem<TPZCompElH1<pzshape::TPZShapePoint > >(mesh,gel);
	//return new TPZCompElWithMem< TPZIntelGen< pzshape::TPZShapePoint > >(mesh,gel,index);
}

TPZCompEl * TPZElastoPlasticAnalysis::CreatePrismElWithMem(TPZGeoEl *gel, TPZCompMesh &mesh, int64_t &index)
{
	return new TPZCompElWithMem<TPZCompElH1<pzshape::TPZShapePrism > >(mesh,gel);
	//return new TPZCompElWithMem< TPZIntelGen< pzshape::TPZShapePrism > >(mesh,gel,index);
}

TPZCompEl * TPZElastoPlasticAnalysis::CreatePyramElWithMem(TPZGeoEl *gel, TPZCompMesh &mesh, int64_t &index)
{
	return new TPZCompElWithMem<TPZCompElH1<pzshape::TPZShapePiram > >(mesh,gel);
	//return new TPZCompElWithMem< TPZIntelGen< pzshape::TPZShapePiram > >(mesh,gel,index);
}

TPZCompEl * TPZElastoPlasticAnalysis::CreateQuadElWithMem(TPZGeoEl *gel, TPZCompMesh &mesh, int64_t &index)
{
//	return new TPZCompElWithMem< TPZIntelGenPlus<TPZIntelGen< pzshape::TPZShapeQuad > > >(mesh,gel,index);
	return new TPZCompElWithMem<TPZCompElH1<pzshape::TPZShapeQuad > >(mesh,gel);
	//return new TPZCompElWithMem< TPZIntelGen< pzshape::TPZShapeQuad > > (mesh,gel,index);
}

TPZCompEl * TPZElastoPlasticAnalysis::CreateTetraElWithMem(TPZGeoEl *gel, TPZCompMesh &mesh, int64_t &index)
{
	return new TPZCompElWithMem<TPZCompElH1<pzshape::TPZShapeTetra > >(mesh,gel);
	//return new TPZCompElWithMem< TPZIntelGen< pzshape::TPZShapeTetra > >(mesh,gel,index);
}

TPZCompEl * TPZElastoPlasticAnalysis::CreateTriangElWithMem(TPZGeoEl *gel, TPZCompMesh &mesh, int64_t &index)
{
	return new TPZCompElWithMem<TPZCompElH1<pzshape::TPZShapeTriang > >(mesh,gel);
	//return new TPZCompElWithMem< TPZIntelGen< pzshape::TPZShapeTriang > >(mesh,gel,index);
}


void TPZElastoPlasticAnalysis::IdentifyEquationsToZero()
{
    fEquationstoZero.clear();
    int64_t nel = fCompMesh->NElements();
    for (int64_t iel=0; iel<nel; iel++) {
        TPZCompEl *cel = fCompMesh->ElementVec()[iel];
        if (!cel) {
            continue;
        }
        TPZMaterial *mat = cel->Material();
        if (!mat) {
            continue;
        }
        int matid = mat->Id();
        if (fMaterialIds.find(matid) == fMaterialIds.end()) {
            continue;
        }
        std::pair<std::multimap<int, int>::iterator,std::multimap<int, int>::iterator> ret;
        ret = fMaterialIds.equal_range(matid);
        std::multimap<int, int>::iterator it;
        for (it=ret.first; it != ret.second; it++)
        {
            int direction = it->second;
            int64_t nc = cel->NConnects();
            for (int64_t ic=0; ic<nc; ic++) {
                TPZConnect &c = cel->Connect(ic);
                int64_t seqnum = c.SequenceNumber();
                int64_t pos = fCompMesh->Block().Position(seqnum);
                int blsize = fCompMesh->Block().Size(seqnum);
                for (int64_t i=pos+direction; i<pos+blsize; i+=2) {
                    fEquationstoZero.insert(i);
                }
            }
        }
    }
#ifdef PZ_LOG
    {
        if(EPAnalysisLogger.isDebugEnabled())
        {
            std::stringstream sout;
            sout << "Equations to zero ";
            std::set<int64_t>::iterator it;
            for (it=fEquationstoZero.begin(); it!= fEquationstoZero.end(); it++) {
                sout << *it << " ";
            }
            LOGPZ_DEBUG(EPAnalysisLogger, sout.str())
        }
    }
#endif
}

/// return the vector of active equation indices
void TPZElastoPlasticAnalysis::GetActiveEquations(TPZVec<int64_t> &activeEquations)
{
    int64_t neq = fCompMesh->NEquations();
    TPZVec<int> equationflag(neq,1);
    typedef std::set<int64_t>::iterator setit;
    for (setit it = fEquationstoZero.begin(); it != fEquationstoZero.end(); it++) {
        equationflag[*it] = 0;
    }
    activeEquations.resize(neq-fEquationstoZero.size());
    int64_t count = 0;
    for (int64_t i=0; i<neq; i++) {
        if (equationflag[i]==1) {
            activeEquations[count++] = i;
        }
    }
}

void  TPZElastoPlasticAnalysis::LoadSolution ( TPZFMatrix<STATE> & loadsol )
{
    fSolution = loadsol;
    LoadSolution();
}


