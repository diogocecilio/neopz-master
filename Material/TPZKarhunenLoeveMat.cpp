#include "TPZKarhunenLoeveMat.h"


TPZKarhunenLoeveMat::TPZKarhunenLoeveMat ( int id,REAL Lx,REAL Ly,REAL Lz, int dim , int type,int expansionorder )
{
  fId=id;
  fLx=Lx;
  fLy=Ly;
  fLz=Lz;
  fDim=dim;
  ftype=type;
  fExpasionOrder=expansionorder;
  SetId ( id );
}


TPZKarhunenLoeveMat::TPZKarhunenLoeveMat ( int id ) :fId ( id )
{
  SetId ( id );
}

TPZKarhunenLoeveMat::TPZKarhunenLoeveMat()
{

}

TPZKarhunenLoeveMat::TPZKarhunenLoeveMat ( const TPZKarhunenLoeveMat &mat ) :fId ( mat.fId ),fLx ( mat.fLx ),fLy ( mat.fLy ),fLz(mat.fLz),ftype ( mat.ftype )
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
  return TPZKarhunenLoeveMat::NSolutionVariables ( var );
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

REAL TPZKarhunenLoeveMat::AutocorrelationFunc ( TPZManVector<REAL,3>  x1, TPZManVector<REAL,3>  x2 )
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
    case 3:
      val = exp ( -fabs ( xx1 - xx2 ) / ( fLx )- fabs ( yy1 - yy2 ) / ( fLy ) );
      break;
	case 4:
      val = exp ( -dist/ fLx );
      break;
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

  TPZVec<STATE> &Sol =data.sol[0];
        switch ( var )
        {
        case EVEC:
          
          //Sol =data.sol[0];
          Solout[0] = Sol[0];
          break;
        case EVEC1:
          Sol =data.sol[1];
          Solout[0] = Sol[0];
          break;
                  case EVEC2:
          Sol =data.sol[2];
          Solout[0] = Sol[0];
          break;
                  case EVEC3:
          Sol =data.sol[3];
          Solout[0] = Sol[0];
          break;
                  case EVEC4:
          Sol =data.sol[4];
          Solout[0] = Sol[0];
          break;
                  case EVEC5:
          Sol =data.sol[5];
          Solout[0] = Sol[0];
          break;
        case EVECSQR:
          Sol =data.sol[0];
          Solout[0] = Sol[0]*Sol[0];
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
