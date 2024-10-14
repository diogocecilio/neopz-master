//
//
#include "TPZRandomFieldAnalysis.h"
#include "TPZEigenSolver.h"
#include "TPZKrylovEigenSolver.h"
#include "TPZSpStructMatrix.h"
#include "pzysmp.h"
#include "pzsysmp.h"
#include "TPZSkylineNSymStructMatrix.h"
#include "TPZMatGeneralisedEigenVal.h"
#include "TPZMaterial.h"
#include "pzcmesh.h"
#include "pzlog.h"
#include <random>
#ifdef PZ_LOG
static TPZLogger logger("pz.analysis");
static TPZLogger loggerError("pz.analysis.error");
#endif

TPZRandomFieldAnalysis::TPZRandomFieldAnalysis(): TPZRegisterClassId(&TPZRandomFieldAnalysis::ClassId),TPZAnalysis()
{

}

TPZRandomFieldAnalysis::TPZRandomFieldAnalysis(TPZCompMesh *mesh,
                                   bool mustOptimizeBandwidth, std::ostream &out): TPZAnalysis(mesh, mustOptimizeBandwidth,out)
{

}

template<class TVar>
TPZEigenSolver<TVar> &TPZRandomFieldAnalysis::EigenSolver()
{
    const auto tmp = dynamic_cast<TPZEigenSolver<TVar>*>(fSolver);
    PZError<<" Aborting...\n";
    DebugStop();
    return *tmp;
}

void TPZRandomFieldAnalysis::SetSolver(const TPZSolver &solver)
{

    PZError<<" Aborting...\n";
    DebugStop();
}


void TPZRandomFieldAnalysis::Assemble()
{

  if(!fCompMesh){
    std::stringstream sout;
    sout<<__PRETTY_FUNCTION__;
    sout<<"\nNo computational mesh found!\n";
    return;
  }

  TPZFMatrix<REAL> Bmat,Cmat;

  AssembleB(Bmat);

  AssembleC(Cmat);

  TPZVec<CSTATE> eigenvalues;
  TPZFMatrix<CSTATE> eigenvectors;

  TPZLapackEigenSolver<REAL> defaultSolver;
 //  defaultSolver.SetNEigenpairs(-60);

   defaultSolver.SolveGeneralisedEigenProblem(Cmat,Bmat,eigenvalues,eigenvectors);

   for(int i=0;i<eigenvalues.size();i++)cout<< eigenvalues[i] <<endl;

   LoadRealSol(eigenvectors);

   TPZFMatrix<REAL> normalizedeigenvectors = NormalizeSolution();

  int ncols =normalizedeigenvectors.Cols();
  int nrows = normalizedeigenvectors.Rows();
  cout << "ncols " << ncols << endl;
  cout << "nrows " << nrows << endl;

 // normalizedeigenvectors.Print("eigenvectors");

  TPZFMatrix<REAL> solut(nrows,ncols);

  for(int ivalue=0;ivalue<nrows;ivalue++)
  {
    //cout << eigenvalues[ivalue].real()  <<endl;
    for(int isol=0;isol<ncols;isol++)
    {
      solut ( ivalue,isol ) =sqrt (eigenvalues[isol].real() ) *normalizedeigenvectors( ivalue, isol );
    }
  }
  fCompMesh->LoadSolution(solut);


}


void TPZRandomFieldAnalysis::Solve()
{

      TPZVec<TPZFMatrix<REAL>> fEigenVectors;

      TPZFMatrix<REAL> fEigenValues;

        std::ofstream outfull ( "outinfofulltime-p2-h4.txt" );
        chrono::steady_clock fulltime;
        auto startfull = fulltime.now();
        std::ofstream out ( "outinfo-p2-h4.txt" );
        TPZFMatrix<REAL> invB,invBC,B,C;
        MatrixXd eigenInvBC,eigenInvB,eigenC,eigenB;
        cout << "Number of Equations =   " << fCompMesh->NEquations() << std::endl;
        cout << "Assembling C  " << std::endl;
        chrono::steady_clock sc;
        auto start = sc.now();

        AssembleC ( C );

        //int nrows = C.Rows();
        //int ncols = C.Rows();

        auto end = sc.now();
        auto time_span = static_cast<chrono::duration<double>> ( end - start );
        cout << "| total time taken to assemble C =  " << time_span.count() << std::endl;

        cout << " Assembling B  " << std::endl;
        start = sc.now();

        AssembleB ( B );

        end = sc.now();
        time_span = static_cast<chrono::duration<double>> ( end - start );
        cout << "| total time taken to assemble B =  " << time_span.count() << std::endl;

        ToEigen ( B,eigenB );
        ToEigen ( C,eigenC );

        eigenInvB = eigenB.inverse();
        eigenInvBC = eigenInvB*eigenC;

        MatrixXd val,vec;

        //GeneralizedEigenSolver<MatrixXd> ces;

        //EigenSolver<MatrixXd> ces;
        ComplexEigenSolver<MatrixXd> ces;

        cout << " Computing Eigensystem  " << std::endl;
        start = sc.now();
        ces.compute ( eigenInvBC );
        //ces.compute(eigenB, eigenC);
        end = sc.now();
        time_span = static_cast<chrono::duration<double>> ( end - start );
        cout << "| total time taken to compute the Eigensystem =  " << time_span.count() << std::endl;
        ces.info();

        auto endfull = fulltime.now();
        time_span = static_cast<chrono::duration<double>> ( endfull - startfull );
        cout << "| total time  =  " << time_span.count() << std::endl;

        int ncols = ces.eigenvectors().cols();
        if(!fNValues||fNValues>ncols)
        {
          fNValues = ncols;
        }
        int nrows = ces.eigenvectors().rows();
        val.resize ( fNValues,1 );
        vec.resize ( nrows,fNValues);

        for ( int icol=0; icol< fNValues; icol++ ) {
                val ( icol,0 ) =ces.eigenvalues() [nrows-icol-1].real();
                for ( int irow=0; irow<nrows; irow++ ) {
                        vec ( irow,icol ) =ces.eigenvectors() ( irow,ncols-icol-1 ).real();
                }
        }

        bool print =true;
        if ( print==true ) {
                cout << "Eigenvalues "<<endl;
                cout << val << endl;
                //cout << "Eigenvectors "<<endl;
                //cout << vec << endl;
                // std::cout << " A "<<endl;
                // std::cout << eigenInvBC <<std::endl;
                //std::cout << " A- V * D * V^(-1) = " << "\n" << (vec * val.asDiagonal() * vec.inverse()) << endl;
        }

        cout << "Start to integrate the solution over the domain..  " << std::endl;
        start = sc.now();
        TPZFMatrix<REAL> vecpz ( nrows,fNValues );
        fEigenVectors.Resize ( fNValues );
        VectorXd intphisqr ( fNValues );
        for ( int i=0; i<fNValues; i++ ) fEigenVectors[i].Resize ( nrows,1 );
        for ( int icol=0; icol< fNValues; icol++ ) {
                for ( int irow=0; irow<nrows; irow++ ) {
                        fEigenVectors[icol] ( irow,0 ) =vec.col ( icol ) ( irow );
                }

                //fEigenVectors[icol].Print(std::cout);
                LoadSolution ( fEigenVectors[icol] );
                int varid=6;
                REAL integral = IntegrateSolution ( varid );
                intphisqr ( icol ) =integral;
                fEigenVectors[icol]*=1./sqrt ( integral );
                vec.col ( icol ) *=1./sqrt ( integral );
        }
        end = sc.now();
        time_span = static_cast<chrono::duration<double>> ( end - start );
        cout << "| total time taken to integrate the solution over the domain =  " << time_span.count() << std::endl;

        //cout << intphisqr << endl;
        for ( int icol=0; icol< fNValues; icol++ ) {
                for ( int irow=0; irow<nrows; irow++ ) {
                        vecpz ( irow,icol ) =sqrt ( val ( icol ) ) *vec.col ( icol ) ( irow );
                }
        }
        LoadSolution ( vecpz );


}


void TPZRandomFieldAnalysis::AssembleC (TPZFMatrix<double> &C)
{
  long iel, jel;
  long nelem =  fCompMesh->NElements();
  TPZManVector<TPZManVector<int>> meshtopology;

  TPZElementMatrixT<REAL> temp (  fCompMesh, TPZElementMatrix::EK );
  TPZAdmChunkVector<TPZCompEl *> &elementvec =  fCompMesh->ElementVec();
  int sz =  fCompMesh->NEquations();
  C.Resize ( sz, sz );

  long count = 0;
  for ( iel = 0; iel < nelem; iel++ )
    {
      //std::cout << "\n iel " << iel << std::endl;
      TPZCompEl *el = elementvec[iel];
      TPZInterpolationSpace *intel = dynamic_cast<TPZInterpolationSpace*>(el);
      if(!intel)
      {
        DebugStop();
      }
      for ( int jel = 0; jel < nelem; jel++ )
        {
          //std::cout << "\n jel " << jel << std::endl;
          TPZCompEl *elj = elementvec[jel];
          TPZElementMatrixT<REAL> ce (  fCompMesh, TPZElementMatrix::EK );
          intel->CalcStiffC (elj, ce );
          //ce.fMat.Print(std::cout);
          int nshape = ce.fMat.Rows();
          TPZManVector<long> SourceIndexIEL, DestinationIndexIEL, SourceIndexJEL, DestinationIndexJEL;
          GetDestIndex ( iel, nshape, SourceIndexIEL, DestinationIndexIEL );
          GetDestIndex ( jel, nshape, SourceIndexJEL, DestinationIndexJEL );

          for ( int irow = 0; irow < DestinationIndexIEL.size(); irow++ )
            {
              for ( int icol = 0; icol < DestinationIndexJEL.size(); icol++ )
                {
                  //REAL val=  C->GetVal( DestinationIndexIEL[irow], DestinationIndexJEL[icol] );
                  //val+=ce.fMat ( SourceIndexIEL[irow], SourceIndexJEL[icol] );
                  //C->PutVal( DestinationIndexIEL[irow], DestinationIndexJEL[icol],val);
                  C ( DestinationIndexIEL[irow], DestinationIndexJEL[icol] ) += ce.fMat ( SourceIndexIEL[irow], SourceIndexJEL[icol] );
                }
            }

        }//jel

    }//iel

  //C.Print ( std::cout );
}

void TPZRandomFieldAnalysis::AssembleB(TPZFMatrix<double> &B)
{
  long iel;
  long nelem = fCompMesh->NElements();
  TPZManVector<TPZManVector<int>> meshtopology;


  TPZElementMatrixT<REAL> temp ( fCompMesh, TPZElementMatrix::EK );
  TPZAdmChunkVector<TPZCompEl *> &elementvec = fCompMesh->ElementVec();
  int64_t sz = fCompMesh->NEquations();
  B.Resize(sz,sz);

  for ( iel = 0; iel < nelem; iel++ )
    {
      TPZCompEl *el = elementvec[iel];
      TPZElementMatrixT<REAL> be ( fCompMesh, TPZElementMatrix::EK );
      TPZInterpolationSpace *intel = dynamic_cast<TPZInterpolationSpace*>(el);
      if(!intel)
      {
        DebugStop();
      }
      intel->CalcStiffB ( be );
      int nshape = be.fMat.Rows();
      TPZManVector<long> SourceIndexIEL, DestinationIndexIEL;
      GetDestIndex ( iel, nshape, SourceIndexIEL, DestinationIndexIEL );
      for ( int irow = 0; irow < DestinationIndexIEL.size(); irow++ )
          {
            for ( int icol = 0; icol < DestinationIndexIEL.size(); icol++ )
              {
                  //REAL val=  B->GetVal( DestinationIndexIEL[irow], DestinationIndexIEL[icol]  );
                 // val+=be.fMat ( SourceIndexIEL[irow], SourceIndexIEL[icol] );
                 // B->PutVal( DestinationIndexIEL[irow], DestinationIndexIEL[icol],val);
                  B ( DestinationIndexIEL[irow], DestinationIndexIEL[icol] ) += be.fMat ( SourceIndexIEL[irow], SourceIndexIEL[icol] );
              }
          }
    }//iel

  //B.Print ( std::cout );
}


void TPZRandomFieldAnalysis::GetDestIndex ( int iel, int nshape, TPZManVector<long> &fSourceIndex, TPZManVector<long> &fDestinationIndex )
{
  fSourceIndex.Resize ( nshape );
  fDestinationIndex.Resize ( nshape );
  long destindex = 0L;
  long fullmatindex = 0L;
  const int numnod = fCompMesh->ElementVec() [iel]->NConnects();
  for ( int in = 0; in < numnod; in++ )
    {
      const long npindex = fCompMesh->ElementVec() [iel]->ConnectIndex ( in );
      TPZConnect &np = fCompMesh->ConnectVec() [npindex];
      long blocknumber = np.SequenceNumber();
      long firsteq = fCompMesh->Block().Position ( blocknumber );
      int ndf = fCompMesh->Block().Size ( blocknumber );
      if ( np.HasDependency() || np.IsCondensed() )
        {
          fullmatindex += ndf;
          continue;
        }//for (np)
      for ( int idf = 0; idf < ndf; idf++ )
        {
          fSourceIndex[destindex] = fullmatindex++;
          fDestinationIndex[destindex++] = firsteq + idf;
        }//for idf
    }//for in

}

TPZFMatrix<REAL> TPZRandomFieldAnalysis::NormalizeSolution()
{
        cout << "Start to integrate the solution over the domain..  " << std::endl;

        TPZFMatrix<REAL> solcopy(fSolution);
        int nrows= fSolution.Rows();
        int ncols= fSolution.Cols();
        int varid=6;
        REAL integral=0.;
        for (int irow=0;irow<nrows;irow++ )
        {
          for (int icol=0; icol< ncols; icol++)
          {
                integral = IntegrateSolution ( varid );
                REAL temp =solcopy.GetVal(irow,icol);
                temp*=1./sqrt ( integral );
                solcopy.PutVal(irow,icol,temp);
          }
        }
        return solcopy;
}

REAL TPZRandomFieldAnalysis::IntegrateSolution ( int varid )
{
        int nvar=1;
        int nels = fCompMesh->NElements(), iel;
        REAL sum=0.;
        TPZVec<REAL> vecsol ( nvar );
        for ( iel = 0; iel < nels; iel++ ) {
                TPZCompEl * cel = fCompMesh->Element ( iel );
                vecsol = cel->IntegrateSolution ( varid );
                sum+=vecsol[0];
        }
        return  sum ;
}
void TPZRandomFieldAnalysis::LoadRealSol(TPZFMatrix<CSTATE> sol)
{
  TPZFMatrix<REAL> solut(sol.Rows(),sol.Cols());
  for(int ivalue=0;ivalue<sol.Rows();ivalue++)
  {
    //cout << sol[ivalue].real()  <<endl;
    for(int isol=0;isol<sol.Cols();isol++)
    {
      solut ( ivalue,isol ) =sol( ivalue, isol ).real();
    }
  }
  LoadSolution(solut);
}
void TPZRandomFieldAnalysis::ManageFieldCretion()
{
  int nfields = fMeanvec.size();
  fFields.resize(nfields);
  if(!nfields)DebugStop();
  for(int ifield=0;ifield<nfields;ifield++)
  {
    fFields[ifield] = GenerateNonGaussinRandomField ( fMeanvec[ifield],fCovvec[ifield],fSamples );
  }

}
TPZFMatrix<REAL>  TPZRandomFieldAnalysis::GenerateNonGaussinRandomField (REAL mean, REAL cov,int samples )
{
        TPZVec<REAL> mean2 ( 2 );
        TPZVec<REAL> cov2 ( 2 );
        TPZVec<string> file2 ( 2 );


        TPZFMatrix<REAL>  PHIt,PHI=fSolution,soltemp=fSolution;

        int M = PHI.Cols();

        TPZFMatrix<REAL> THETA ( M,samples );

        std::normal_distribution<REAL> distribution ( 0., 1. );
        for ( int n = 0; n < samples; n++ )
        {
          for ( int iexp = 0; iexp < M; iexp++ )
          {
            std::random_device rd{};
            std::mt19937 generator{ rd() };
            REAL xic = distribution ( generator );
            THETA ( iexp,n ) = xic;
          }
        }

        TPZFMatrix<REAL> hhat;
        PHI.Multiply ( THETA, hhat );

        REAL sdev = cov* mean;
        REAL xi = sqrt ( log ( 1 + pow ( ( sdev / mean),2 ) ) );
        REAL lambda = log ( mean) - xi * xi / 2.;

        for ( int i = 0; i < hhat.Rows(); i++ )
        {
          for ( int j = 0; j < hhat.Cols(); j++ )
          {
              hhat ( i,j ) = exp ( lambda + xi * hhat ( i,j ) );
          }
        }

        //PrintField(hhat,filename);
        //LoadSolution(hhat);
        return hhat;
}

void TPZRandomFieldAnalysis::PrintField(TPZFMatrix<REAL> matfield,string filename)
{
          std::ofstream print ( filename );
          int rows = matfield.Rows();
          int cols = matfield.Cols();
          for(int irow=0;irow<rows;irow++)
          {
            for(int icol=0;icol<cols;icol++)
            {
              print << matfield(irow,icol) << " ";
            }
            print << "\n";
          }

}

void TPZRandomFieldAnalysis::FromEigen ( MatrixXd eigenmat, TPZFMatrix<REAL>  &pzmat )
    {
		int rows = eigenmat.rows();
		int cols = eigenmat.cols();
		pzmat.Resize(rows,cols);
		for(int irow=0;irow<rows;irow++)
		{
			for(int icol=0;icol<cols;icol++)
			{
				pzmat(irow,icol)=eigenmat(irow,icol);
			}
		}

    }

void TPZRandomFieldAnalysis::ToEigen ( TPZFMatrix<REAL>  pzmat,MatrixXd &eigenmat )
    {
		TPZFMatrix<REAL> intpz(pzmat);
		int rows = pzmat.Rows();
		int cols = pzmat.Cols();
		eigenmat.resize(rows,cols);
		for(int irow=0;irow<rows;irow++)
		{
			for(int icol=0;icol<cols;icol++)
			{
				eigenmat(irow,icol)=intpz(irow,icol);
			}
		}

    }

  int TPZRandomFieldAnalysis::ClassId() const
{
  return Hash("TPZRandomFieldAnalysis") ^
    TPZAnalysis::ClassId() << 1;
}

void TPZRandomFieldAnalysis::Write(TPZStream &buf, int withclassid) const
{
  //TPZAnalysis::Write(buf,withclassid);
  fSolution.Write(buf,withclassid);
  for(int ifield=0;ifield< fFields.size(); ifield++)
  {
    fFields[ifield].Write(buf,withclassid);
  }


}

void TPZRandomFieldAnalysis::Read(TPZStream &buf, void *context)
{
  //TPZAnalysis::Read(buf,context);
  fSolution.Read(buf,context);
  fFields.resize(2);
  for(int ifield=0;ifield< fFields.size(); ifield++)
  {
    fFields[ifield].Read(buf,context);
  }

}
