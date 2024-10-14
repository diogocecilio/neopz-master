#include "TPZKLStrMatrix.h"
#include "pzfstrmatrix.h"
#include "pzfmatrix.h"
#include "pzcmesh.h"
#include "pzsubcmesh.h"
#include <sstream>
#include "pzlog.h"
#include "pzfmatrix.h"
#include "pzmatrix.h"
#include "tpzklinterpolationspace.h"

#include "pzfstrmatrix.h"
#include "pzfmatrix.h"
#include "pzsubcmesh.h"
#include <sstream>
#include "pzlog.h"

using namespace std;


using namespace std;
template<class TVar, class TPar>
TPZMatrix<TVar> * TPZKLStrMatrix<TVar,TPar>::Create(){
	int64_t neq = this->fEquationFilter.NActiveEquations();

	return new TPZFMatrix<TVar>(neq,neq,0.);
}

template<class TVar, class TPar>
TPZStructMatrix * TPZKLStrMatrix<TVar,TPar>::Clone(){
    return new TPZKLStrMatrix(*this);
}


template<class TVar, class TPar>
int TPZKLStrMatrix<TVar,TPar>::ClassId() const{
    return Hash("TPZKLStrMatrix") ^
        TPZStructMatrixT<TVar>::ClassId() << 1 ^
        TPar::ClassId() << 2;
}

template<class TVar, class TPar>
void TPZKLStrMatrix<TVar,TPar>::Read(TPZStream& buf, void* context){
    TPZStructMatrix::Read(buf,context);
    TPar::Read(buf,context);
}

template<class TVar, class TPar>
void TPZKLStrMatrix<TVar,TPar>::Write(TPZStream& buf, int withclassid) const{
    TPZStructMatrix::Write(buf,withclassid);
    TPar::Write(buf,withclassid);
}

template<class TVar, class TPar>
void TPZKLStrMatrix<TVar,TPar>::AssembleC (TPZFMatrix<REAL> &C)
{
  long iel, jel;
  long nelem =  this->fMesh->NElements();
  TPZManVector<TPZManVector<int>> meshtopology;

  TPZElementMatrixT<REAL> temp (  this->fMesh, TPZElementMatrix::EK );
  TPZAdmChunkVector<TPZCompEl *> &elementvec =  this->fMesh->ElementVec();
  int sz =  this->fMesh->NEquations();
  C.Resize ( sz, sz );

  long count = 0;
  for ( iel = 0; iel < nelem; iel++ )
    {
      std::cout << "\n iel " << iel << std::endl;
      TPZCompEl *el = elementvec[iel];
      TPZKLInterpolationSpace *intel = dynamic_cast<TPZKLInterpolationSpace*>(el);
      if(!intel)
      {
        DebugStop();
      }
      for ( int jel = 0; jel < nelem; jel++ )
        {
          std::cout << "\n jel " << jel << std::endl;
          TPZCompEl *elj = elementvec[jel];
          TPZElementMatrixT<REAL> ce (  this->fMesh, TPZElementMatrix::EK );
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
                  C ( DestinationIndexIEL[irow], DestinationIndexJEL[icol] ) += ce.fMat ( SourceIndexIEL[irow], SourceIndexJEL[icol] );
                }
            }

        }//jel

    }//iel

  C.Print ( std::cout );
}

template<class TVar, class TPar>
void TPZKLStrMatrix<TVar,TPar>::AssembleB(TPZFMatrix<REAL> &B)
{
  long iel;
  long nelem = this->fMesh->NElements();
  TPZManVector<TPZManVector<int>> meshtopology;

  TPZElementMatrixT<REAL> temp ( this->fMesh, TPZElementMatrix::EK );
  TPZAdmChunkVector<TPZCompEl *> &elementvec = this->fMesh->ElementVec();
  int sz = this->fMesh->NEquations();
  B.Resize ( sz, sz );

  for ( iel = 0; iel < nelem; iel++ )
    {
      TPZCompEl *el = elementvec[iel];
      TPZElementMatrixT<REAL> be ( this->fMesh, TPZElementMatrix::EK );
      TPZKLInterpolationSpace *intel = dynamic_cast<TPZKLInterpolationSpace*>(el);
      intel->CalcStiffB ( be );
      int nshape = be.fMat.Rows();
      TPZManVector<long> SourceIndexIEL, DestinationIndexIEL;
      GetDestIndex ( iel, nshape, SourceIndexIEL, DestinationIndexIEL );
      for ( int irow = 0; irow < DestinationIndexIEL.size(); irow++ )
          {
            for ( int icol = 0; icol < DestinationIndexIEL.size(); icol++ )
              {
                B ( DestinationIndexIEL[irow], DestinationIndexIEL[icol] ) += be.fMat ( SourceIndexIEL[irow], SourceIndexIEL[icol] );
              }
          }
    }//iel

 // C.Print ( std::cout );
}
template<class TVar, class TPar>
void TPZKLStrMatrix<TVar,TPar>::GetDestIndex ( int iel, int nshape, TPZManVector<long> &fSourceIndex, TPZManVector<long> &fDestinationIndex )
{
  fSourceIndex.Resize ( nshape );
  fDestinationIndex.Resize ( nshape );
  long destindex = 0L;
  long fullmatindex = 0L;
  const int numnod = this->fMesh->ElementVec() [iel]->NConnects();
  for ( int in = 0; in < numnod; in++ )
    {
      const long npindex = this->fMesh->ElementVec() [iel]->ConnectIndex ( in );
      TPZConnect &np = this->fMesh->ConnectVec() [npindex];
      long blocknumber = np.SequenceNumber();
      long firsteq = this->fMesh->Block().Position ( blocknumber );
      int ndf = this->fMesh->Block().Size ( blocknumber );
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

#include "pzstrmatrixot.h"
#include "pzstrmatrixflowtbb.h"

template class TPZKLStrMatrix<STATE,TPZStructMatrixOR<STATE>>;
template class TPZKLStrMatrix<STATE,TPZStructMatrixOT<STATE>>;
template class TPZKLStrMatrix<STATE,TPZStructMatrixTBBFlow<STATE>>;

template class TPZKLStrMatrix<CSTATE,TPZStructMatrixOR<CSTATE>>;
template class TPZKLStrMatrix<CSTATE,TPZStructMatrixOT<CSTATE>>;
template class TPZKLStrMatrix<CSTATE,TPZStructMatrixTBBFlow<CSTATE>>;

