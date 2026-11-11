#include "TPZKarhunenLoeveMat.h"


TPZKarhunenLoeveMat::TPZKarhunenLoeveMat ( int id,REAL Lx,REAL Ly,REAL Lz, int dim , int type )
{
  fId=id;
  fLx=Lx;
  fLy=Ly;
  fLz=Lz;
  fDim=dim;
  ftype=type;

  SetId ( id );
}


// the constructors without the kernel parameters leave them with defined values (before, they were not initialized):
// dimension 2 and type 0 (no kernel, AutocorrelationFunc stops)
TPZKarhunenLoeveMat::TPZKarhunenLoeveMat ( int id ) :fId ( id ), fLx ( 1. ), fLy ( 1. ), fLz ( 1. ), ftype ( 0 ), fDim ( 2 )
{
  SetId ( id );
}

TPZKarhunenLoeveMat::TPZKarhunenLoeveMat() : fId ( 0 ), fLx ( 1. ), fLy ( 1. ), fLz ( 1. ), ftype ( 0 ), fDim ( 2 )
{

}

// fDim is copied too (before, the copy had it not initialized)
TPZKarhunenLoeveMat::TPZKarhunenLoeveMat ( const TPZKarhunenLoeveMat &mat ) :fId ( mat.fId ),fLx ( mat.fLx ),fLy ( mat.fLy ),fLz(mat.fLz),ftype ( mat.ftype ), fDim ( mat.fDim )
{
  SetId ( mat.fId );
}

TPZKarhunenLoeveMat::~TPZKarhunenLoeveMat()
{


}



int TPZKarhunenLoeveMat::VariableIndex(const std::string &name) const
{
  if ( !strcmp ( "vec",		name.c_str() ) )  return TPZKarhunenLoeveMat::EVEC;
  if ( !strcmp ( "vec1",		name.c_str() ) )  return TPZKarhunenLoeveMat::EVEC1;
  if ( !strcmp ( "vec2",		name.c_str() ) )  return TPZKarhunenLoeveMat::EVEC2;
  if ( !strcmp ( "vec3",		name.c_str() ) )  return TPZKarhunenLoeveMat::EVEC3;
  if ( !strcmp ( "vec4",		name.c_str() ) )  return TPZKarhunenLoeveMat::EVEC4;
  if ( !strcmp ( "vec5",		name.c_str() ) )  return TPZKarhunenLoeveMat::EVEC5;
  if ( !strcmp ( "vecsqr",		name.c_str() ) )  return TPZKarhunenLoeveMat::EVECSQR;
  PZError << "TPZMatElastoPlastic::VariableIndex Error\n";
//DebugStop();
  return -1;
}


int TPZKarhunenLoeveMat::NSolutionVariables(int var) const
{
  if ( var == TPZKarhunenLoeveMat::EVEC )		 return 1;
  if ( var == TPZKarhunenLoeveMat::EVEC1 )		 return 1;
  if ( var == TPZKarhunenLoeveMat::EVEC2 )		 return 1;
  if ( var == TPZKarhunenLoeveMat::EVEC3 )		 return 1;
  if ( var == TPZKarhunenLoeveMat::EVEC4 )		 return 1;
  if ( var == TPZKarhunenLoeveMat::EVEC5 )		 return 1;
  if ( var == TPZKarhunenLoeveMat::EVECSQR )		 return 1;
  // the base class reports the unknown variable (before, this method called itself: infinite recursion)
  return TBase::NSolutionVariables ( var );
}



void TPZKarhunenLoeveMat::ContributeB ( TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ek )
{
  TPZFMatrix<REAL>  &phi = data.phi;
  TPZFMatrix<REAL>  phit;
  phi.Transpose ( &phit );
  phi.Multiply ( phit, ek );

  ek *= weight;
}

void TPZKarhunenLoeveMat::ContributeC ( TPZMaterialData &data1,TPZMaterialData &data2, REAL w1,REAL w2, TPZFMatrix<STATE> &ek )
{

  TPZFMatrix<REAL>  &phi1 = data1.phi;
  TPZFMatrix<REAL>  &phi2 = data2.phi;
  TPZFMatrix<REAL>  phi1t,phi2t;

  phi2.Transpose ( &phi2t );

  REAL CXX = AutocorrelationFunc ( data1.x,data2.x );

  phi1.Multiply ( phi2t, ek );
  ek *= CXX*w1*w2;
}

void TPZKarhunenLoeveMat::ContributeC ( const TPZKLPointData &data1, const TPZKLPointData &data2, TPZFMatrix<STATE> &ek ) const
{
  const int64_t n1 = data1.fPhi.Rows(), np1 = data1.fPhi.Cols();
  const int64_t n2 = data2.fPhi.Rows(), np2 = data2.fPhi.Cols();
  if ( data1.fW.size() != np1 || data2.fW.size() != np2 || data1.fX.Cols() != np1 || data2.fX.Cols() != np2 )
    {
      PZError << __PRETTY_FUNCTION__ << " inconsistent point data\n";
      DebugStop();
    }
  ek.Redim ( n1, n2 );
  if ( n1 == 0 || n2 == 0 || np1 == 0 || np2 == 0 ) return;
  // sum over the points q of the second element first: g_j(x_p) = sum_q G(x_p,x_q) w2_q phi2_j(x_q),
  // then ek(i,j) += w1_p phi1_i(x_p) g_j(x_p). The same result as phi1 phi2^T G w1 w2 summed point pair by point pair,
  // with O(np1 np2 n2 + np1 n1 n2) operations
  TPZManVector<REAL,3> x1 ( 3 ), x2 ( 3 );
  TPZManVector<STATE,64> g ( n2 );
  const REAL *phi1 = &data1.fPhi.g ( 0, 0 ); // column major: phi1 + p * n1 is the column of the point p
  const REAL *phi2 = &data2.fPhi.g ( 0, 0 );
  STATE *ekp = &ek ( 0, 0 );
  for ( int64_t p = 0; p < np1; p++ )
    {
      for ( int k = 0; k < 3; k++ ) x1[k] = data1.fX.GetVal ( k, p );
      for ( int64_t j = 0; j < n2; j++ ) g[j] = 0.;
      for ( int64_t q = 0; q < np2; q++ )
        {
          for ( int k = 0; k < 3; k++ ) x2[k] = data2.fX.GetVal ( k, q );
          const REAL gw = AutocorrelationFunc ( x1, x2 ) * data2.fW[q];
          const REAL *phi2q = phi2 + q * n2;
          for ( int64_t j = 0; j < n2; j++ ) g[j] += gw * phi2q[j];
        }
      const REAL w1 = data1.fW[p];
      const REAL *phi1p = phi1 + p * n1;
      for ( int64_t j = 0; j < n2; j++ )
        {
          const STATE gj = w1 * g[j];
          STATE *col = ekp + j * n1;
          for ( int64_t i = 0; i < n1; i++ ) col[i] += phi1p[i] * gj;
        }
    }
}

REAL TPZKarhunenLoeveMat::AutocorrelationFunc ( const TPZVec<REAL> &x1, const TPZVec<REAL> &x2 ) const
{
  REAL val = 0, xx1, xx2, yy1, yy2,zz1,zz2, dist;

  xx1 = x1[0];
  yy1 = x1[1];
  zz1 = x1[2];
  
  xx2 = x2[0];
  yy2 = x2[1];
  zz2 = x2[2];

  dist = sqrt ( pow ( xx1 - xx2, 2 ) + pow ( yy1 - yy2, 2 )+pow ( zz1 - zz2, 2 ) );

  switch ( ftype )
    {
    case 1:
      val =  exp ( -fabs ( xx1 - xx2 ) * fabs ( xx1 - xx2 ) / ( fLx * fLx ) - fabs ( yy1 - yy2 ) * fabs ( yy1 - yy2 ) / ( fLy * fLy ) );
      break;
    case 2:
      val = exp ( - ( ( dist / fLx ) * ( dist / fLx ) ) );
      break; // before, case 2 fell into case 3 and returned the kernel of type 3
    case 3:
      val = exp ( -fabs ( xx1 - xx2 ) / ( fLx )- fabs ( yy1 - yy2 ) / ( fLy ) );
      break;
	case 4:
      val = exp ( -dist/ fLx );
      break;
    default:
      PZError << __PRETTY_FUNCTION__ << " autocorrelation type " << ftype << " not defined (1 to 4)\n";
      DebugStop();
    }
  if ( val>1 )
    {
      std::cout <<"The AutocorrelationFunc cant be larger than 1" <<std::endl;
      DebugStop();
    }


  return val;
}

void TPZKarhunenLoeveMat::Contribute(const TPZMaterialDataT<STATE> &data,
                                 REAL weight,
                                 TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef)
{
  DebugStop();
}



void TPZKarhunenLoeveMat::ContributeBC(const TPZMaterialDataT<STATE> &data,
                                   REAL weight,
                                   TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef,
                                   TPZBndCondT<STATE> &bc)
{
  DebugStop();
}

void TPZKarhunenLoeveMat::Solution(const TPZMaterialDataT<STATE> &data,
                               int var, TPZVec<STATE> &Solout)
{

  // before, Sol was a reference to data.sol[0] and "Sol = data.sol[k]" overwrote data.sol[0] with the load case k
        switch ( var )
        {
        case EVEC:
          Solout[0] = data.sol[0][0];
          break;
        case EVEC1:
          Solout[0] = data.sol[1][0];
          break;
        case EVEC2:
          Solout[0] = data.sol[2][0];
          break;
        case EVEC3:
          Solout[0] = data.sol[3][0];
          break;
        case EVEC4:
          Solout[0] = data.sol[4][0];
          break;
        case EVEC5:
          Solout[0] = data.sol[5][0];
          break;
        case EVECSQR:
          Solout[0] = data.sol[0][0]*data.sol[0][0];
          break;
        }
        
 }
 void TPZKarhunenLoeveMat::Errors(const TPZMaterialDataT<STATE> &data,
                             TPZVec<REAL> &values) {
   std::cout << "Not implemented" <<std::endl;
   DebugStop();
 }
int TPZKarhunenLoeveMat::ClassId() const {
    return Hash("TPZKarhunenLoeveMat") ^ TBase::ClassId() << 1;
}

void TPZKarhunenLoeveMat::FillDataRequirements(TPZMaterialData &data) const
{
    data.fNeedsSol = true;
    data.fNeedsNormal = false;

}

void TPZKarhunenLoeveMat::FillBoundaryConditionDataRequirements(int type,
                                                           TPZMaterialData &data) const
{
    data.fNeedsSol = false;
    data.fNeedsNormal = false;
    if (type == 4 || type == 5 || type == 6) {
        data.fNeedsNormal = true;
    }
}
int TPZKarhunenLoeveMat::NStateVariables() const {
	return 1;
}
