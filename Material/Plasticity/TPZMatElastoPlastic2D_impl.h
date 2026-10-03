//
//  TPZMatElastoPlastic2D_impl.h
//  pz
//
//  Created by Omar Durán on 9/8/18.
//

#include "TPZMatElastoPlastic2D.h"
#include "TPZMatElastoPlastic_impl.h"

#include "TPZBndCond.h"

#ifndef WIN32
#include <fenv.h>//NAN DETECTOR
#endif

#ifdef PZ_LOG
#include "pzlog.h"
static TPZLogger elastoplastic2dLogger("material.pzElastoPlastic2D");
#endif


template <class T, class TMEM>
TPZMatElastoPlastic2D<T,TMEM>::TPZMatElastoPlastic2D() : TPZMatElastoPlastic<T,TMEM>()
{
    fPlaneStrain = true;
}

template <class T, class TMEM>
TPZMatElastoPlastic2D<T,TMEM>::TPZMatElastoPlastic2D(int id , int PlaneStrainOrPlaneStress) : TPZMatElastoPlastic<T,TMEM>(id)
{
    fPlaneStrain = PlaneStrainOrPlaneStress;
}

template <class T, class TMEM>
TPZMatElastoPlastic2D<T,TMEM>::TPZMatElastoPlastic2D(const TPZMatElastoPlastic2D<T,TMEM> &mat) : TPZMatElastoPlastic<T,TMEM>(mat), fPlaneStrain(mat.fPlaneStrain)
{
#ifdef PZ_LOG
    if(elastoplastic2dLogger.isDebugEnabled())
    {
        std::stringstream sout;
        sout << ">>> TPZMatElastoPlastic2D<T,TMEM>() copy constructor called ***";
        LOGPZ_INFO(elastoplastic2dLogger,sout.str().c_str());
    }
#endif
}

template <class T, class TMEM>
TPZMatElastoPlastic2D<T,TMEM>::~TPZMatElastoPlastic2D()
{
    
}

template <class T, class TMEM>
void TPZMatElastoPlastic2D<T,TMEM>::ApplyDeltaStrain(const TPZMaterialDataT<STATE> & data, TPZFMatrix<REAL> & DeltaStrain,TPZFMatrix<REAL> & Stress)
{
    
    
#ifdef PZDEBUG
    if (DeltaStrain.Rows() != 6) {
        DebugStop();
    }
#endif
    
    if (!fPlaneStrain) //
    {//
        DebugStop();//PlaneStress
    }
    
    TPZMatElastoPlastic<T,TMEM>::ApplyDeltaStrain(data,DeltaStrain,Stress);//
    
}


template <class T, class TMEM>
void TPZMatElastoPlastic2D<T,TMEM>::ApplyDeltaStrainComputeDep(const TPZMaterialDataT<STATE> & data, TPZFMatrix<REAL> & DeltaStrain,TPZFMatrix<REAL> & Stress, TPZFMatrix<REAL> & Dep)
{
#ifdef PZDEBUG
    if (DeltaStrain.Rows() != 6) {
        DebugStop();
    }
#endif
    if (!fPlaneStrain) //
    {//
        DebugStop();//PlaneStress
    }
    TPZMatElastoPlastic<T,TMEM>::ApplyDeltaStrainComputeDep(data,DeltaStrain,Stress,Dep);//
}

template <class T, class TMEM>
void TPZMatElastoPlastic2D<T, TMEM>::Contribute(const TPZMaterialDataT<STATE> &data, REAL weight, TPZFMatrix<REAL> &ek, TPZFMatrix<REAL> &ef) {
    
    const TPZFMatrix<REAL> &dphi = data.dphix;
    const TPZFMatrix<REAL> &phi = data.phi;
    const TPZFMatrix<REAL> &axes = data.axes;
    TPZFMatrix<REAL> dphiXY, axesT;
    axes.Transpose(&axesT);
    axesT.Multiply(dphi, dphiXY);
    
    const int phr = phi.Rows();
    
    TPZFNMatrix<4> Deriv(2, 2);
    TPZFNMatrix<36> Dep(6, 6,0.0);
    TPZFNMatrix<6> DeltaStrain(6, 1);
    TPZFNMatrix<6> Stress(6, 1);
    int ptindex = data.intGlobPtIndex;
    
    if (TPZMatWithMem<TMEM>::fUpdateMem && data.sol.size() > 0) {
        // Loop over the solutions if update memory is true
        TPZSolVec<STATE> locsol(data.sol);
        TPZGradSolVec<STATE> locdsol(data.dsol);
        int numsol = locsol.size();
        
        for (int is = 0; is < numsol; is++) {
            data.sol[0] = locsol[is];
            data.dsol[0] = locdsol[is];
            
            this->ComputeDeltaStrainVector(data, DeltaStrain);
            this->ApplyDeltaStrainComputeDep(data, DeltaStrain, Stress, Dep);
        }
    } else {
        this->ComputeDeltaStrainVector(data, DeltaStrain);
        this->ApplyDeltaStrainComputeDep(data, DeltaStrain, Stress, Dep);
    }
    
    ptindex = 0;
    int nstate = NStateVariables();
    REAL val;
    TPZManVector<STATE,3> ForceLoc(this->m_force);
    //TPZManVector<STATE, 2> ForceLoc(nstate,0.0);
    if (this->fForcingFunction) {
        this->fForcingFunction(data.x, ForceLoc);
    }

    // ek += w Bᵀ D B e ef += w (Nᵀ f - Bᵀ σ), com B em Voigt de engenharia (ε_xx, ε_yy, γ_xy) — mesma conta de
    // BuildBu/Dep2D, escrita com laços para não alocar matrizes temporárias em cada ponto de integração
    const STATE D[3][3] = {{Dep(_XX_, _XX_), Dep(_XX_, _YY_), Dep(_XX_, _XY_)},
                           {Dep(_YY_, _XX_), Dep(_YY_, _YY_), Dep(_YY_, _XY_)},
                           {Dep(_XY_, _XX_), Dep(_XY_, _YY_), Dep(_XY_, _XY_)}};
    const STATE sxx = Stress(_XX_, 0), syy = Stress(_YY_, 0), sxy = Stress(_XY_, 0);
    for (int in = 0; in < phr; in++) {
        const STATE dxi = dphiXY(0, in), dyi = dphiXY(1, in);
        // linhas de Bᵀ para os graus de liberdade (in, u) e (in, v): Bu = (dx, 0, dy), Bv = (0, dy, dx)
        const STATE Bu[3] = {dxi, 0., dyi}, Bv[3] = {0., dyi, dxi};
        STATE DBu[3], DBv[3];  // (Bᵀ D) das duas linhas
        for (int k = 0; k < 3; k++) {
            DBu[k] = Bu[0] * D[0][k] + Bu[1] * D[1][k] + Bu[2] * D[2][k];
            DBv[k] = Bv[0] * D[0][k] + Bv[1] * D[1][k] + Bv[2] * D[2][k];
        }
        for (int jn = 0; jn < phr; jn++) {
            const STATE dxj = dphiXY(0, jn), dyj = dphiXY(1, jn);
            // colunas de B: Bu_j = (dxj, 0, dyj), Bv_j = (0, dyj, dxj)
            ek(nstate * in + 0, nstate * jn + 0) += weight * (DBu[0] * dxj + DBu[2] * dyj);
            ek(nstate * in + 0, nstate * jn + 1) += weight * (DBu[1] * dyj + DBu[2] * dxj);
            ek(nstate * in + 1, nstate * jn + 0) += weight * (DBv[0] * dxj + DBv[2] * dyj);
            ek(nstate * in + 1, nstate * jn + 1) += weight * (DBv[1] * dyj + DBv[2] * dxj);
        }
        ef(nstate * in + 0, 0) += weight * (ForceLoc[0] * phi(in, 0) - (dxi * sxx + dyi * sxy));
        ef(nstate * in + 1, 0) += weight * (ForceLoc[1] * phi(in, 0) - (dyi * syy + dxi * sxy));
    }
}

template <class T, class TMEM>
void TPZMatElastoPlastic2D<T, TMEM>::Contribute(const TPZMaterialDataT<STATE> &data, REAL weight, TPZFMatrix<REAL> &ef) {
    const TPZFMatrix<REAL> &dphi = data.dphix;
    const TPZFMatrix<REAL> &phi = data.phi;
    const TPZFMatrix<REAL> &axes = data.axes;
    TPZFNMatrix<9, REAL> axesT;
    TPZFNMatrix<50, REAL> dphiXY;
    
    axes.Transpose(&axesT);
    axesT.Multiply(dphi, dphiXY);
    
    const int phr = phi.Rows();
    
    //TPZFNMatrix<36> Deriv(6, 6);
    TPZFNMatrix<6> DeltaStrain(6, 1);
    TPZFNMatrix<6> Stress(6, 1);
    int ptindex = data.intGlobPtIndex;
    
    if (TPZMatWithMem<TMEM>::fUpdateMem && data.sol.size() > 0) {
        
        TPZSolVec<STATE> locsol(data.sol);
        TPZGradSolVec<STATE> locdsol(data.dsol);
        int numsol = locsol.size();
        
        for (int is = 0; is < numsol; is++) {
            data.sol[0] = locsol[is];
            data.dsol[0] = locdsol[is];
            
            this->ComputeDeltaStrainVector(data, DeltaStrain);
            this->ApplyDeltaStrain(data, DeltaStrain, Stress);
        }
    } else {
        this->ComputeDeltaStrainVector(data, DeltaStrain);
        this->ApplyDeltaStrain(data, DeltaStrain, Stress);
    }
#ifdef MACOS
    feclearexcept(FE_ALL_EXCEPT);
    if (fetestexcept(/*FE_DIVBYZERO*/ FE_ALL_EXCEPT)) {
        std::cout << "division by zero reported\n";
        DebugStop();
    }
#endif
    
#ifdef PZ_LOG
    if (elastoplastic2dLogger.isDebugEnabled()) {
        std::stringstream sout;
        sout << ">>> TPZMatElastoPlastic<T,TMEM>::Contribute ***";
        sout << "\nIntegration Local Point index = " << data.intGlobPtIndex;
        sout << "\nIntegration Global Point index = " << data.intGlobPtIndex;
        sout << "\ndata.axes = " << data.axes;
        sout << "\nStress \n";
        sout << Stress(_XX_, 0) << "\t" << Stress(_YY_, 0) << "\t" << Stress(_XY_, 0) << "\n";
        sout << "\nDELTA STRAIN \n";
        sout << DeltaStrain(0, 0) << "\t" << DeltaStrain(1, 0) << "\t" << DeltaStrain(2, 0) << "\n";
        sout << "data.phi" << data.phi<< std::endl;
        TPZMatWithMem<TMEM>::MemItem(ptindex).Print(sout);
        
        LOGPZ_DEBUG(elastoplastic2dLogger, sout.str().c_str());
    }
#endif
    int nstate = NStateVariables();
    REAL val;
    
    //TPZManVector<STATE, 2> ForceLoc(nstate,0.0);
    TPZManVector<STATE,3> ForceLoc(this->m_force);
    if (this->fForcingFunction) {
        this->fForcingFunction(data.x, ForceLoc);
    }


    int in;
	for(in = 0; in < phr; in++)
	{
        val  =((ForceLoc[0]) * phi(in,0));
        val -= Stress(_XX_,0) * dphiXY(0,in);
        val -= Stress(_XY_,0) * dphiXY(1,in);
        ef(in*nstate+0,0) += weight * val;

        val  = ((ForceLoc[1]) * phi(in,0));
        val -= Stress(_XY_,0) * dphiXY(0,in);
        val -= Stress(_YY_,0) * dphiXY(1,in);
        ef(in*nstate+1,0) += weight * val;
    }
    
}


template <class T, class TMEM>
void TPZMatElastoPlastic2D<T,TMEM>::FillBoundaryConditionDataRequirements(int type,TPZMaterialData &data) const
{
    
    data.fNeedsSol = true;
    if (type == 4 || type ==5 || type == 6) {
        data.fNeedsNormal = true;
    }
    else {
        data.fNeedsNormal = false;
    }
}

template <class T, class TMEM>
void TPZMatElastoPlastic2D<T, TMEM>::ContributeBC(const TPZMaterialDataT<STATE> &data, REAL weight, TPZFMatrix<REAL> &ek, TPZFMatrix<REAL> &ef, TPZBndCondT<STATE> &bc) {
    int nstate = NStateVariables();
    const REAL BIGNUMBER = TPZMaterial::fBigNumber;

    auto &bc_with_memory = dynamic_cast<TPZMatWithMem<TMEM> &>(bc);  // referência: antes copiava toda a memória do contorno em cada ponto

    /// Accepting  solution on bc data.
    int gp_index = data.intGlobPtIndex;
    const TPZFMatrix<REAL> &phi = data.phi;
    TPZManVector<STATE,3> delta_u    = data.sol[0];
    TPZManVector<STATE,3> u_n(nstate,0.0);
    TPZManVector<STATE,3> u(bc_with_memory.MemItem(gp_index).m_u);
    for (int i = 0; i < nstate; i++) {
        u_n[i] = delta_u[i] + u[i];
    }
    
    if(TPZMatWithMem<TMEM>::fUpdateMem)
    {
        bc_with_memory.MemItem(gp_index).m_u = u_n;
    }
    
    const int phr = phi.Rows();
    int in, jn, idf, jdf;
    const auto &v2 = bc.Val2();
    
    const TPZFMatrix<REAL> &v1 = bc.Val1();
    switch (bc.Type()) {
        case 0: // Dirichlet condition
            for (in = 0; in < phr; in++) {
                ef(nstate * in + 0, 0) += BIGNUMBER * (v2[0] - u_n[0]) * phi(in, 0) * weight;
                ef(nstate * in + 1, 0) += BIGNUMBER * (v2[1] - u_n[1]) * phi(in, 0) * weight;
                
                for (jn = 0; jn < phr; jn++) {
                    ek(nstate * in + 0, nstate * jn + 0) += BIGNUMBER * phi(in, 0) * phi(jn, 0) * weight;
                    ek(nstate * in + 1, nstate * jn + 1) += BIGNUMBER * phi(in, 0) * phi(jn, 0) * weight;
                    
                }//jn
            }//in
            break;
            
        case 1: // Neumann condition
            for (in = 0; in < phi.Rows(); in++) {

              ef(nstate*in+0,0) += (v2[0] * phi(in,0) * weight);
              ef(nstate*in+1,0) += (v2[1] * phi(in,0) * weight);
            }
        break;
        case 2: // Mixed condition (mola): mesma forma do ContributeBC só de resíduo; phi é (nshape x 1)
        {
            TPZFNMatrix<2, STATE> res(2, 1, 0.);
            for (int i = 0; i < 2; i++) {
                for (int j = 0; j < 2; j++) {
                    res(i, 0) += bc.Val1()(i, j) * u_n[j];
                }
            }

            for (in = 0; in < phr; in++) {
                ef(nstate * in + 0, 0) += (v2[0] - res(0, 0)) * phi(in, 0) * weight;
                ef(nstate * in + 1, 0) += (v2[1] - res(1, 0)) * phi(in, 0) * weight;
                for (jn = 0; jn < phr; jn++) {
                    for (idf = 0; idf < 2; idf++) {
                        for (jdf = 0; jdf < 2; jdf++) {
                            ek(nstate * in + idf, nstate * jn + jdf) += bc.Val1()(idf, jdf) * phi(in, 0) * phi(jn, 0) * weight;
                        }
                    }
                }
            }// este caso pode reproduzir o caso 0 quando o deslocamento
        }
        break;


        case 3: // Directional Null Dirichlet - displacement is set to null in the non-null vector component direction
            for (in = 0; in < phr; in++) {
                ef(nstate * in + 0, 0) += BIGNUMBER * (0. - u_n[0]) * v2[0] * phi(in, 0) * weight;
                ef(nstate * in + 1, 0) += BIGNUMBER * (0. - u_n[1]) * v2[1] * phi(in, 0) * weight;
                for (jn = 0; jn < phr; jn++) {
                    ek(nstate * in + 0, nstate * jn + 0) += BIGNUMBER * phi(in, 0) * phi(jn, 0) * weight * v2[0];
                    ek(nstate * in + 1, nstate * jn + 1) += BIGNUMBER * phi(in, 0) * phi(jn, 0) * weight * v2[1];
                }//jn
            }//in
            break;
            
        case 4: // stressField Neumann condition
        {
            // cópia local: antes escrevia em bc.Val2() através da referência const (alterava o contorno)
            TPZManVector<STATE,2> v2loc(2, 0.);
            v2loc[0] = v1(0, 0) * data.normal[0] + v1(0, 1) * data.normal[1];
            v2loc[1] = v1(1, 0) * data.normal[0] + v1(1, 1) * data.normal[1];
            // The normal vector points towards the neighbor. The negative sign is there to
            // reflect the outward normal vector.
            for (in = 0; in < phi.Rows(); in++) {
                ef(nstate * in + 0, 0) += v2loc[0] * phi(in, 0) * weight;
                ef(nstate * in + 1, 0) += v2loc[1] * phi(in, 0) * weight;
            }
        }
            break;
        case 5://PRESSAO DEVE SER POSTA NA POSICAO 0 DO VETOR v2
        {
            TPZFNMatrix<2, STATE> res(2, 1, 0.);
            for (int i = 0; i < 2; i++) {
                for (int j = 0; j < 2; j++) {
                    res(i, 0) += data.normal[i] * bc.Val1()(i, j) * u_n[j] * data.normal[j];
                }
            }
            for (in = 0; in < phi.Rows(); in++) {
                ef(nstate * in + 0, 0) += (v2[0] * data.normal[0] - res(0, 0)) * phi(in, 0) * weight;
                ef(nstate * in + 1, 0) += (v2[0] * data.normal[1] - res(1, 0)) * phi(in, 0) * weight;
                for (jn = 0; jn < phi.Rows(); jn++) {
                    for (idf = 0; idf < 2; idf++) {
                        for (jdf = 0; jdf < 2; jdf++) {
                            ek(nstate * in + idf, nstate * jn + jdf) += bc.Val1()(idf, jdf) * data.normal[idf] * data.normal[jdf] * phi(in, 0) * phi(jn, 0) * weight;
                            // BUG FALTA COLOCAR VAL2
                            // DebugStop();
                        }
                    }
                }
                
            }
        }
            break;
            
        case 6:
        {
            
            REAL v[1];
            v[0] = bc.Val2()[0];    //    Tn normal component of normal traction (T)
            
            REAL tn = v[0];
            TPZManVector<REAL,3> n = data.normal;
            ///    Neumann condition for each state variable
            ///    Elasticity Equation
            for(in = 0 ; in < phi.Rows(); in++)
            {
                ///   Normal Tension Components on neumman boundary
                ef(nstate*in+0,0)    += weight * tn * n[0] * phi(in,0);        //    Tnx
                ef(nstate*in+1,0)    += weight * tn * n[1] * phi(in,0);        //    Tny
            }
            
        }
            break;
            
        case 7 : {//colocar 1 em val1(1,1)=1 para impor deslocamento val2(1)=uprescribedy
            
            REAL v_null[2];
            v_null[0] = bc.Val1()(0, 0);
            v_null[1] = bc.Val1()(1, 1);
            
            for (in = 0; in < phr; in++) {
                ef(nstate * in + 0, 0) += BIGNUMBER * (v2[0] - u_n[0]) * v_null[0] * phi(in, 0) * weight;
                ef(nstate * in + 1, 0) += BIGNUMBER * (v2[1] - u_n[1]) * v_null[1] * phi(in, 0) * weight;
                for (jn = 0; jn < phr; jn++) {
                    ek(nstate * in + 0, nstate * jn + 0) += BIGNUMBER * phi(in, 0) * phi(jn, 0) * weight * v_null[0];
                    ek(nstate * in + 1, nstate * jn + 1) += BIGNUMBER * phi(in, 0) * phi(jn, 0) * weight * v_null[1];
                }//jn
            }//in
            break;
        }

        default:
#ifdef PZ_LOG
            if (elastoplastic2dLogger.isDebugEnabled()) {
                std::stringstream sout;
                sout << "<<< TPZMatElastoPlastic2D<T,TMEM>::ContributeBC *** WRONG BOUNDARY CONDITION TYPE = " << bc.Type();
                LOGPZ_ERROR(elastoplastic2dLogger, sout.str().c_str());
            }
#endif
            PZError << "TPZMatElastoPlastic2D::ContributeBC error - Wrong boundary condition type" << std::endl;
    }
}

template <class T, class TMEM>
void TPZMatElastoPlastic2D<T, TMEM>::ContributeBC(const TPZMaterialDataT<STATE> &data, REAL weight, TPZFMatrix<STATE> &ef, TPZBndCondT<STATE> &bc) {
    TPZFMatrix<STATE> ek_fake(ef.Rows(),ef.Rows());
    int nstate = NStateVariables();
    const REAL BIGNUMBER = TPZMaterial::fBigNumber;
    
    auto * bc_with_memory = dynamic_cast<TPZMatWithMem<TMEM> *>(&bc);
    if(!bc_with_memory){
        PZError << "TPZMatElastoPlastic2D::ContributeBC error - Wrong boundary objected: expected to be TPZMatWithMem<TMEM>" << std::endl;
        DebugStop();
    }
    
    /// Accepting  solution on bc data.
    int gp_index = data.intGlobPtIndex;
    const TPZFMatrix<REAL> &phi = data.phi;
    TPZManVector<STATE,3> delta_u    = data.sol[0];
    TPZManVector<STATE,3> u_n(nstate,0.0);
    TPZManVector<STATE,3> u(bc_with_memory->MemItem(gp_index).m_u);
    for (int i = 0; i < nstate; i++) {
        u_n[i] = delta_u[i] + u[i];
    }
    
    if(TPZMatWithMem<TMEM>::fUpdateMem)
    {
        bc_with_memory->MemItem(gp_index).m_u = u_n;
    }
    
    const int phr = phi.Rows();
    int in, jn, idf, jdf;
    const auto &v2 = bc.Val2();
    
    const TPZFMatrix<REAL> &v1 = bc.Val1();
    switch (bc.Type()) {
        case 0: // Dirichlet condition
            for (in = 0; in < phr; in++) {
                ef(nstate * in + 0, 0) += BIGNUMBER * (v2[0] - u_n[0]) * phi(in, 0) * weight;
                ef(nstate * in + 1, 0) += BIGNUMBER * (v2[1] - u_n[1]) * phi(in, 0) * weight;
                
            }//in
            break;
            
        case 1: // Neumann condition
            for (in = 0; in < phi.Rows(); in++) {
                ef(nstate * in + 0, 0) += v2[0] * phi(in, 0) * weight;
                ef(nstate * in + 1, 0) += v2[1] * phi(in, 0) * weight;
            }
            break;
            
        case 2: // Mixed condition
        {
            TPZFNMatrix<2, STATE> res(2, 1, 0.);
            for (int i = 0; i < 2; i++) {
                for (int j = 0; j < 2; j++) {
                    res(i, 0) += bc.Val1()(i, j) * u_n[j];
                }
            }
            
            for (in = 0; in < phi.Rows(); in++) {
                ef(nstate * in + 0, 0) += (v2[0] - res(0, 0)) * phi(in, 0) * weight;
                ef(nstate * in + 1, 0) += (v2[1] - res(1, 0)) * phi(in, 0) * weight;
            }//in
        }
            break;
            
        case 3: // Directional Null Dirichlet - displacement is set to null in the non-null vector component direction
            for (in = 0; in < phr; in++) {
                ef(nstate * in + 0, 0) += BIGNUMBER * (0. - u_n[0]) * v2[0] * phi(in, 0) * weight;
                ef(nstate * in + 1, 0) += BIGNUMBER * (0. - u_n[1]) * v2[1] * phi(in, 0) * weight;
            }//in
            break;
            
        case 4: // stressField Neumann condition
        {
            // cópia local: antes escrevia em bc.Val2() através da referência const (alterava o contorno)
            TPZManVector<STATE,2> v2loc(2, 0.);
            v2loc[0] = v1(0, 0) * data.normal[0] + v1(0, 1) * data.normal[1];
            v2loc[1] = v1(1, 0) * data.normal[0] + v1(1, 1) * data.normal[1];
            // The normal vector points towards the neighbor. The negative sign is there to
            // reflect the outward normal vector.
            for (in = 0; in < phi.Rows(); in++) {
                ef(nstate * in + 0, 0) += v2loc[0] * phi(in, 0) * weight;
                ef(nstate * in + 1, 0) += v2loc[1] * phi(in, 0) * weight;
            }
        }
            break;
        case 5://PRESSAO DEVE SER POSTA NA POSICAO 0 DO VETOR v2
        {
            TPZFNMatrix<2, STATE> res(2, 1, 0.);
            for (int i = 0; i < 2; i++) {
                for (int j = 0; j < 2; j++) {
                    res(i, 0) += data.normal[i] * bc.Val1()(i, j) * u_n[j] * data.normal[j];
                }
            }
            for (in = 0; in < phi.Rows(); in++) {
                ef(nstate * in + 0, 0) += (v2[0] * data.normal[0] - res(0, 0)) * phi(in, 0) * weight;
                ef(nstate * in + 1, 0) += (v2[0] * data.normal[1] - res(1, 0)) * phi(in, 0) * weight;
            }
        }
            break;
            
        case 6:
        {
            
            REAL v[1];
            v[0] = bc.Val2()[0];    //    Tn normal component of normal traction (T)
            
            REAL tn = v[0];
            TPZManVector<REAL,3> n = data.normal;
            ///    Neumann condition for each state variable
            ///    Elasticity Equation
            for(in = 0 ; in < phi.Rows(); in++)
            {
                ///   Normal Tension Components on neumman boundary
                ef(nstate*in+0,0)    += weight * tn * n[0] * phi(in,0);        //    Tnx
                ef(nstate*in+1,0)    += weight * tn * n[1] * phi(in,0);        //    Tny
            }
            
        }
            break;
            
        case 7 : {
            
            REAL v_null[2];
            v_null[0] = bc.Val1()(0, 0);
            v_null[1] = bc.Val1()(1, 1);
            
            for (in = 0; in < phr; in++) {
                ef(nstate * in + 0, 0) += BIGNUMBER * (v2[0] - u_n[0]) * v_null[0] * phi(in, 0) * weight;
                ef(nstate * in + 1, 0) += BIGNUMBER * (v2[1] - u_n[1]) * v_null[1] * phi(in, 0) * weight;
            }//in
            break;
            
        }
            
        default:
#ifdef PZ_LOG
            if (elastoplastic2dLogger.isDebugEnabled()) {
                std::stringstream sout;
                sout << "<<< TPZMatElastoPlastic2D<T,TMEM>::ContributeBC *** WRONG BOUNDARY CONDITION TYPE = " << bc.Type();
                LOGPZ_ERROR(elastoplastic2dLogger, sout.str().c_str());
            }
#endif
            PZError << "TPZMatElastoPlastic2D::ContributeBC error - Wrong boundary condition type" << std::endl;
    }
}



template <class T, class TMEM>
void TPZMatElastoPlastic2D<T,TMEM>::Solution(const TPZMaterialDataT<STATE> &data, int var, TPZVec<REAL> &Solout)
{
    
    TPZMaterialDataT<STATE> datalocal(data);
    datalocal.sol[0].Resize(3,0.);
    datalocal.dsol[0].Resize(3,3);
    datalocal.dsol[0](2,0) = 0.;
    datalocal.dsol[0](2,1) = 0.;
    datalocal.dsol[0](2,2) = 0.;
    datalocal.dsol[0](0,2) = 0.;
    datalocal.dsol[0](1,2) = 0.;

    TPZMatElastoPlastic<T,TMEM>::Solution(datalocal,var,Solout);






}



template <class T, class TMEM>
void TPZMatElastoPlastic2D<T, TMEM>::ComputeDeltaStrainVector(const TPZMaterialDataT<STATE> & data, TPZFMatrix<REAL> &DeltaStrain) {
    TPZFNMatrix<9> DSolXYZ(3, 3, 0.);
    data.axes.Multiply(data.dsol[0], DSolXYZ, 1/*transpose*/);
    if (DeltaStrain.Rows() != 6) {
        DebugStop();
    }
    
    DeltaStrain(_XX_, 0) = DSolXYZ(0, 0);
    DeltaStrain(_YY_, 0) = DSolXYZ(1, 1);
    DeltaStrain(_XY_, 0) = (DSolXYZ(1, 0) + DSolXYZ(0, 1));
    DeltaStrain(_XZ_, 0) = 0.;
    DeltaStrain(_YZ_, 0) = 0.;
    DeltaStrain(_ZZ_, 0) = 0.;
    
}


template <class T, class TMEM>
TPZMaterial * TPZMatElastoPlastic2D<T,TMEM>::NewMaterial() const
{
    return new TPZMatElastoPlastic2D<T,TMEM>(*this);
}

#include "TPZSandlerExtended.h"
#include "TPZPlasticStepPV.h"
#include "TPZYCMohrCoulombPV.h"

template <class T, class TMEM>
std::string TPZMatElastoPlastic2D<T,TMEM>::Name() const
{
    return "TPZMatElastoPlastic<T,TMEM>";
}

template <class T, class TMEM>
void TPZMatElastoPlastic2D<T,TMEM>::Write(TPZStream &buf, int withclassid) const{
    TPZMatElastoPlastic<T,TMEM>::Write(buf,withclassid);
    int classid = ClassId();
    buf.Write(&classid);
}

template <class T, class TMEM>
void TPZMatElastoPlastic2D<T,TMEM>::Read(TPZStream &buf, void *context)
{
    TPZMatElastoPlastic<T,TMEM>::Read(buf,context);
    int classid;
    buf.Read(&classid);
    if (classid != ClassId()) {
        DebugStop();
    }
}


template <class T, class TMEM>
void TPZMatElastoPlastic2D<T,TMEM>::Print(std::ostream &out, const int memory) const
{
    out << __PRETTY_FUNCTION__ << std::endl;
    TPZMatElastoPlastic<T,TMEM>::Print(out,memory);
}

template <class T, class TMEM>
void TPZMatElastoPlastic2D<T,TMEM>::Print(std::ostream &out) const
{
    out << __PRETTY_FUNCTION__ << std::endl;
    out << "Plane strain " << fPlaneStrain << std::endl;
    TPZMatElastoPlastic<T,TMEM>::Print(out);
}
