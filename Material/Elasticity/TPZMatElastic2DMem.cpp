#include "TPZMatElastic2DMem.h"
#include "TPZMaterialDataT.h" // definição completa de TPZMaterialDataT

// ----------------- helpers -----------------
template<class TMEM>
void TPZMatElastic2DMem<TMEM>::ElasticD(STATE E, STATE nu, TPZFNMatrix<9,STATE>& D) const
{
    D.Redim(3,3); D.Zero();
    if (fPlaneStress) {
        const STATE c = E/(1.0 - nu*nu);
        D(0,0)= c;     D(0,1)= c*nu;  D(0,2)= 0.;
        D(1,0)= c*nu;  D(1,1)= c;     D(1,2)= 0.;
        D(2,0)= 0.;    D(2,1)= 0.;    D(2,2)= E/(2.0*(1.0+nu));
    } else { // plane strain
        const STATE lambda = E*nu/((1.0+nu)*(1.0-2.0*nu));
        const STATE mu     = E/(2.0*(1.0+nu));
        D(0,0)= lambda+2*mu; D(0,1)= lambda;       D(0,2)= 0.;
        D(1,0)= lambda;      D(1,1)= lambda+2*mu;  D(1,2)= 0.;
        D(2,0)= 0.;          D(2,1)= 0.;          D(2,2)= mu;
    }
}

template<class TMEM>
STATE TPZMatElastic2DMem<TMEM>::EfromMem(const TPZMaterialDataT<STATE>& data)
{
    STATE E = fEref;
    auto mem = this->GetMemory();
    const long ig = data.intGlobPtIndex;
    if (mem && ig >= 0 && ig < mem->NElements()) {
        E = (*mem)[ig].fE;
    }
    return E;
}

// ----------------- requisitos -----------------
template<class TMEM>
void TPZMatElastic2DMem<TMEM>::FillDataRequirements(TPZMaterialDataT<STATE>& data) const
{
    TBase::FillDataRequirements(data);   // defaults
    data.SetAllRequirements(false);
    data.fNeedsSol     = true;
    //data.fNeedsGradSol = true;   // se você calcula tensões
    data.fNeedsHSize   = true;
}
template<class TMEM>
void TPZMatElastic2DMem<TMEM>::FillBoundaryConditionDataRequirements(
    int, TPZMaterialDataT<STATE>& data) const
{
    TBase::FillDataRequirements(data);   // defaults
    data.SetAllRequirements(false);
    data.fNeedsSol     = true;
    //data.fNeedsGradSol = true;   // se você calcula tensões
    data.fNeedsHSize   = true;
}

// ----------------- montagem -----------------
template<class TMEM>
void TPZMatElastic2DMem<TMEM>::Contribute(const TPZMaterialDataT<STATE>& data,
                                          REAL weight, TPZFMatrix<STATE>& ek,
                                          TPZFMatrix<STATE>& ef)
{
    const int nshape = data.phi.Rows();
    if (nshape==0) return;

    const auto &dphix = data.dphix; // [dim x nshape]
    const STATE E  = EfromMem(data);
    //const STATE E  = 1.;
    const STATE nu = fNu;
   // const STATE nu = 0.;

    TPZFNMatrix<9,STATE> D; ElasticD(E, nu, D);

    auto dot3 = [&](const STATE s[3], const STATE t[3]) -> STATE {
        const STATE Dt0 = D(0,0)*t[0] + D(0,1)*t[1] + D(0,2)*t[2];
        const STATE Dt1 = D(1,0)*t[0] + D(1,1)*t[1] + D(1,2)*t[2];
        const STATE Dt2 = D(2,0)*t[0] + D(2,1)*t[1] + D(2,2)*t[2];
        return s[0]*Dt0 + s[1]*Dt1 + s[2]*Dt2;
    };

    for (int a=0; a<nshape; a++) {
        const STATE dNax = dphix(0,a);
        const STATE dNay = dphix(1,a);
        const STATE su[3] = { dNax, 0.,   dNay };
        const STATE sv[3] = { 0.,   dNay, dNax };

        for (int b=0; b<nshape; b++) {
            const STATE dNbx = dphix(0,b);
            const STATE dNby = dphix(1,b);
            const STATE tu[3] = { dNbx, 0.,   dNby };
            const STATE tv[3] = { 0.,   dNby, dNbx };

            const STATE val_uu = dot3(su,tu) * weight;
            const STATE val_uv = dot3(su,tv) * weight;
            const STATE val_vu = dot3(sv,tu) * weight;
            const STATE val_vv = dot3(sv,tv) * weight;

            ek(2*a+0, 2*b+0) += val_uu;
            ek(2*a+0, 2*b+1) += val_uv;
            ek(2*a+1, 2*b+0) += val_vu;
            ek(2*a+1, 2*b+1) += val_vv;
        }
    }
    (void)ef; // sem força de corpo
}

template<class TMEM>
void TPZMatElastic2DMem<TMEM>::ContributeBC(const TPZMaterialDataT<STATE>& data,
                                            REAL weight, TPZFMatrix<STATE>& ek,
                                            TPZFMatrix<STATE>& ef,
                                            TPZBndCondT<STATE>& bc)
{
    const int nshape = data.phi.Rows();
    const auto &phi = data.phi;
    const int type = bc.Type();
    constexpr STATE big = (STATE)1.e12;
    if (type == 0) { // Dirichlet (penalidade)

        const auto &val2 = bc.Val2(); // vetor [u_Dx, u_Dy]
        for (int i=0;i<nshape;i++){
            for (int j=0;j<nshape;j++){
                const STATE k = big * phi(i,0) * phi(j,0) * weight;
                ek(2*i+0,2*j+0) += k;
                ek(2*i+1,2*j+1) += k;
            }
            ef(2*i+0) += big * phi(i,0) * (val2.size()>0 ? val2[0] : (STATE)0.) * weight;
            ef(2*i+1) += big * phi(i,0) * (val2.size()>1 ? val2[1] : (STATE)0.) * weight;
        }
    } else if (type == 1) { // Neumann: tração
        const auto &val2 = bc.Val2(); // vetor [tx, ty]
        const STATE tx = (val2.size()>0 ? val2[0] : (STATE)0.);
        const STATE ty = (val2.size()>1 ? val2[1] : (STATE)0.);
        for (int i=0;i<nshape;i++){
            ef(2*i+0) += phi(i,0) * tx * weight;
            ef(2*i+1) += phi(i,0) * ty * weight;
        }
    } else if (type == 3) { // Directional Dirichlet
        for(int in = 0 ; in < nshape; in++) {
            const auto &v2 = bc.Val2();
            for (int jn = 0 ; jn < nshape; jn++) {
                ek(2*in+0,2*jn+0) += big * phi(in,0) * phi(jn,0) * weight * v2[0];
                ek(2*in+1,2*jn+1) += big * phi(in,0) * phi(jn,0) * weight * v2[1];
            }//jn
        }//in
    } else {
        // outros tipos: sem contribuição especial
    }
}

// ----------------- pós-processo -----------------
template<class TMEM>
void TPZMatElastic2DMem<TMEM>::Solution(const TPZMaterialDataT<STATE>& data,
                                        int var, TPZVec<STATE>& out)
{
    const STATE E  = EfromMem(data);
    const STATE nu = fNu;
    TPZFNMatrix<9,STATE> D; ElasticD(E, nu, D);

    // u = (ux, uy)
    const STATE ux = (data.sol.size() && data.sol[0].size() > 0) ? data.sol[0][0] : (STATE)0.;
    const STATE uy = (data.sol.size() && data.sol[0].size() > 1) ? data.sol[0][1] : (STATE)0.;

    // grad u
    STATE du_dx=0., du_dy=0., dv_dx=0., dv_dy=0.;
    if (data.dsol.size()){
        du_dx = data.dsol[0](0,0);
        du_dy = data.dsol[0](1,0);
        dv_dx = data.dsol[0](0,1);
        dv_dy = data.dsol[0](1,1);
    }

    const STATE exx = du_dx;
    const STATE eyy = dv_dy;
    const STATE gxy = du_dy + dv_dx;

    switch (var){
        case EVar_Displacement:
            out.Resize(2); out[0]=ux; out[1]=uy; break;
        case EVar_Strain:
            out.Resize(3); out[0]=exx; out[1]=eyy; out[2]=gxy; break;
        case EVar_Stress: {
            out.Resize(3);
            out[0] = D(0,0)*exx + D(0,1)*eyy + D(0,2)*gxy; // σx
            out[1] = D(1,0)*exx + D(1,1)*eyy + D(1,2)*gxy; // σy
            out[2] = D(2,0)*exx + D(2,1)*eyy + D(2,2)*gxy; // τxy
        } break;
        case EVar_E:
            out.Resize(1); out[0]=E; break;
        default:
            out.Fill(0.);
    }
}

// instanciamento explícito
template class TPZMatElastic2DMem<TKLPointMem>;
