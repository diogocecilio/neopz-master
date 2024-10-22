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


    int ndof = Cmat.Rows();
    if(!fNValues||fNValues>ndof)
    {
        fNValues = ndof;
    }

      cout << fNValues << endl;
      MatrixXd eigenInvBC,eigenInvB,eigenC,eigenB;
      ToEigen ( Bmat,eigenB );
      ToEigen ( Cmat,eigenC );

        eigenInvB = eigenB.inverse();
        eigenInvBC = eigenInvB*eigenC;

        //cout << eigenInvBC << endl;


  TPZVec<CSTATE> eigenvaluescomplex;
  TPZFMatrix<CSTATE> eigenvectorscomplex;

  TPZLapackEigenSolver<REAL> defaultSolver;
   defaultSolver.SetNEigenpairs(fNValues);

   defaultSolver.SolveGeneralisedEigenProblem(Cmat,Bmat,eigenvaluescomplex,eigenvectorscomplex);

   fEigenvetors = LoadRealSol(eigenvectorscomplex);
   fEigenvalues = LoadRealSol(eigenvaluescomplex);

    for(int i=0;i<fEigenvalues.size();i++)cout<< fEigenvalues[i] <<endl;
    fEigenvetors.Print("vecs");
    NormalizeSolution();
    fSolutionValVec.Resize(ndof,fNValues);
    for(int idof=0;idof<ndof;idof++)
    {
      for(int iM=0;iM<fNValues;iM++)
      {
        fSolutionValVec(idof,iM)=sqrt(fEigenvalues[iM])*fEigenvetors(iM,idof);
      }
    }
   // fSolutionValVec.Print("s");
    LoadSolution(fSolutionValVec);

//    for(int i=0;i<fEigenvalues.size();i++)cout<< fEigenvalues[i] <<endl;
//    fEigenvetors.Print("vecs");
//     TPZFMatrix<REAL> newmat(fNValues,fNValues,0.);
//     for(int icol=0;icol<fNValues;icol++)
//     {
//
//         TPZFMatrix<REAL> ivec(fNValues,1),ivect,mult;
//         for(int irow=0;irow<fEigenvetors.Rows();irow++)ivec(irow,0)=fEigenvetors(irow,icol);
//         ivec.Transpose(&ivect);
//         ivec.Multiply(ivect,mult);
//         mult*=fEigenvalues[icol];
//         newmat+=mult;
//
//
//     }
    //newmat.Print("REC");


}


void TPZRandomFieldAnalysis::Solve()
{

      TPZVec<TPZFMatrix<REAL>> fEigenVectors;

        TPZFMatrix<REAL> invB,invBC,B,C;
        MatrixXd eigenInvBC,eigenInvB,eigenC,eigenB;
        cout << "Number of Equations =   " << fCompMesh->NEquations() << std::endl;
        cout << "Assembling C  " << std::endl;

        AssembleC ( C );

        cout << " Assembling B  " << std::endl;

        AssembleB ( B );

        ToEigen ( B,eigenB );
        ToEigen ( C,eigenC );

        eigenInvB = eigenB.inverse();
        eigenInvBC = eigenInvB*eigenC;

        //cout << eigenC << std::endl;

        MatrixXd eigenvalues,eigenvectors,reconstructed;

        //GeneralizedEigenSolver<MatrixXd> ces;

        //EigenSolver<MatrixXd> ces;
        ComplexEigenSolver<MatrixXd> ces;

        cout << " Computing Eigensystem  " << std::endl;

        ces.compute ( eigenInvBC );
        //ces.compute(eigenB, eigenC);
        ces.info();

        int ncols = ces.eigenvectors().cols();
        if(!fNValues||fNValues>ncols)
        {
          fNValues = ncols;
        }
        int nrows = ces.eigenvectors().rows();
        eigenvalues.resize ( fNValues,1 );
        eigenvectors.resize ( nrows,fNValues);
        fEigenvetors.Resize(nrows,fNValues);

        fEigenvalues.resize(fNValues);
        for ( int icol=0; icol< fNValues; icol++ ) {
                eigenvalues ( icol,0 ) =ces.eigenvalues() [nrows-icol-1].real();
                fEigenvalues[icol]= eigenvalues ( icol,0 );
                for ( int irow=0; irow<nrows; irow++ ) {
                        eigenvectors ( irow,icol ) =ces.eigenvectors() ( irow,ncols-icol-1 ).real();
                }
        }
    cout << eigenvalues << endl;
    cout << eigenvectors << endl;

//         //bool print =true;
//         if ( print==true ) {
//                // cout << "Eigenvalues "<<endl;
//               //  cout << val << endl;
//                 //cout << "Eigenvectors "<<endl;
//                 //cout << vec << endl;
//                 // std::cout << " A "<<endl;
//                 // std::cout << eigenInvBC <<std::endl;
//                 //std::cout << " A- V * D * V^(-1) = " << "\n" << (vec * val.asDiagonal() * vec.inverse()) << endl;
//         }

        cout << "Start to integrate the solution over the domain..  " << std::endl;
        TPZFMatrix<REAL> vecpz ( nrows,fNValues );
        fEigenVectors.Resize ( fNValues );
        VectorXd intphisqr ( fNValues );
        for ( int i=0; i<fNValues; i++ ) fEigenVectors[i].Resize ( nrows,1 );
        for ( int icol=0; icol< fNValues; icol++ ) {
                for ( int irow=0; irow<nrows; irow++ ) {
                        fEigenVectors[icol] ( irow,0 ) =eigenvectors.col ( icol ) ( irow );
                }

                //fEigenVectors[icol].Print(std::cout);
                LoadSolution ( fEigenVectors[icol] );
                int varid=6;
                REAL integral = IntegrateSolution ( varid );
                intphisqr ( icol ) =integral;
                fEigenVectors[icol]*=1./sqrt ( integral );
                eigenvectors.col ( icol ) *=1./sqrt ( integral );
        }
        FromEigen(eigenvectors,fEigenvetors);

        //cout << intphisqr << endl;
        for ( int icol=0; icol< fNValues; icol++ ) {
                for ( int irow=0; irow<nrows; irow++ ) {
                        vecpz ( irow,icol ) =sqrt ( eigenvalues ( icol ) ) *eigenvectors.col ( icol ) ( irow );
                }
        }

        fSolutionValVec=vecpz;
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

void TPZRandomFieldAnalysis::NormalizeSolution()
{

  for ( int icol=0; icol< fEigenvetors.Cols(); icol++ )
  {
                TPZFMatrix<REAL> ivec(fEigenvetors.Rows(),1);
                for(int irow=0;irow<fEigenvetors.Rows();irow++)ivec(irow,0)=fEigenvetors(irow,icol);
                LoadSolution ( ivec );
                int varid=6;
                cout << "nao"<<endl;
                REAL integral = IntegrateSolution ( varid );
                cout << "integral = "<< integral <<endl;
                for(int irow=0;irow<fEigenvetors.Rows();irow++)fEigenvetors(irow,icol)*=1./sqrt ( integral );
  }
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
TPZFMatrix<REAL>  TPZRandomFieldAnalysis::LoadRealSol(TPZFMatrix<CSTATE> sol)
{
  int ncols = sol.Cols();
  int nrows = sol.Rows();
  TPZFMatrix<REAL> solut(sol.Rows(),sol.Cols());
  for ( int icol=0; icol< fNValues; icol++ )
  {
    for ( int irow=0; irow<nrows; irow++ )
    {
      solut ( irow,icol ) =sol( irow,ncols-icol-1 ).real();


    }


  }

  return solut;
}
TPZVec<REAL>  TPZRandomFieldAnalysis::LoadRealSol(TPZVec<CSTATE> sol)
{

  int nrows = sol.size();
  TPZVec<REAL> solut(sol.size());
  for ( int icol=0; icol< fNValues; icol++ )
  {
                solut [ icol] =sol[nrows-icol-1].real();
  }
return solut;
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
  fSolutionValVec.Write(buf,withclassid);
  //output.Write(&seqsize);
  buf.Write(fEigenvalues);
  fEigenvetors.Write(buf,withclassid);
//   for(int ifield=0;ifield< fFields.size(); ifield++)
//   {
//     fFields[ifield].Write(buf,withclassid);
//   }


}

void TPZRandomFieldAnalysis::Read(TPZStream &buf, void *context)
{
  //TPZAnalysis::Read(buf,context);
  fSolution.Read(buf,context);
  fSolutionValVec.Read(buf,context);
  buf.Read(fEigenvalues);
  fEigenvetors.Read(buf,context);
//   fFields.resize(2);
//   for(int ifield=0;ifield< fFields.size(); ifield++)
//   {
//     fFields[ifield].Read(buf,context);
//   }

}
