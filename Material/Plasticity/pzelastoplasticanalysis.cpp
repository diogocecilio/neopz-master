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
#ifdef PZ_LOG
static TPZLogger EPAnalysisLogger("pz.analysis.elastoplastic");
static TPZLogger loggertest("testing");
#endif

using namespace std;


TPZElastoPlasticAnalysis::TPZElastoPlasticAnalysis() : TPZLinearAnalysis(), fPrecond(NULL), fLineSearch(ELineSearch::Dicotomic) {
	//Mesh()->Solution().Zero(); already performed in the nonlinearanalysis base class
	//fSolution.Zero();
}

TPZElastoPlasticAnalysis::TPZElastoPlasticAnalysis(TPZCompMesh *mesh,std::ostream &out, ELineSearch lsearch) : TPZLinearAnalysis(mesh,false), fPrecond(NULL),fLineSearch(lsearch) {

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


bool TPZElastoPlasticAnalysis::FindRoot(int &iters,REAL &resu,REAL &resf)
{
    //METODO DE NEWTON SIMPLES
    // estado e incremento
    TPZFMatrix<STATE> x(Solution()), dx(Solution());
    x.Zero(); dx.Zero();

    const REAL tol   = 1.e-6;
    const int  n_it  = 10;
    const REAL EPS   = 1.e-30; // evita divisão por zero

    //std::cout << "AssembleResidual.."   <<endl;
    // resíduo inicial
    AssembleResidual();
    REAL normrhs0 = Norm(Rhs());
    if (!std::isfinite(normrhs0)) normrhs0 = 1.0;
    if (normrhs0 < EPS) { iters = 0; return true; } // já está resolvido

    REAL normu0 = Norm(Solution());
    if (normu0<=0) normu0 = 1.0;
    //if (normu0 < EPS) { iters = 0; return true; } // já está resolvido

    REAL normrhs = 1.0;
    REAL normdu  = 0.0;
    resf=10000.;
    for (int i = 1; i <= n_it; ++i) {
        // monta tangente e resíduo na solução atual x
        //std::cout << "Assembling.."   <<endl;
        Assemble();

        //std::cout << "Solving.."   <<endl;
        // resolve Δx
        Solve();
        dx = Solution();

        // atualiza solução candidata
        x += dx;
        LoadSolution(x);

        // *** RECOMPUTA o resíduo para a solução ATUALIZADA ***
        AssembleResidual();



        // métricas
        normdu  = Norm(dx);
        normrhs = Norm(Rhs());

        std::cout << "  [it " << i << "] "
        << "||Δu|| = " << normdu
        << " | ||R|| = " << normrhs
        << " | tol = " << tol << std::endl;

        iters = i;
        //if(normrhs>resf)break;
        resu=normdu;
        resf=normrhs;
        // critério de convergência (resíduo relativo)
        if (normrhs < tol &&normdu<tol) {
            std::cout << "normrhs ="<< normrhs<< " normdu ="<< normdu   <<endl;
            return true;
        }

    }
    std::cout << "não convergiu dentro do limite: normrhs ="<< normrhs<< " normdu ="<< normdu   <<endl;
    int numiter=100;
    std::cout << "Tentando line search com  ="<< normrhs<< " tol ="<< tol << " numiter = "<< numiter<<std::endl;
    bool conv = IterativeProcess(std::cout,tol, numiter, true,false,iters);
    //bool conv=false;
    // não convergiu dentro do limite
    if( conv)
    {
        std::cout << "Convergiu com   ="<< iters<< " iteraçoes."<<std::endl;
        return true;
    }else{
        std::cout << "NAO Convergiu com   ="<< iters<< " iteraçoes."<<std::endl;
        return false;
    }

}



bool TPZElastoPlasticAnalysis::NewtonRaphson(bool verbose)
{

    TPZFMatrix<STATE> x(Solution()), dx(Solution());
    x.Zero(); dx.Zero();

    const REAL tol   = 1.e-6;
    const int  n_it  = 20;
    const REAL EPS   = 1.e-30;

    //std::cout << "AssembleResidual.."   <<endl;
    int iters;
    AssembleResidual();

    STATE r0=0.,r1=0.,r2=0.;
    TPZFMatrix<STATE> mr0,mr1,mr2;
    STATE rate0=0.,rate1=0.,rate2=0.;
    REAL normrhs = 1.0;
    REAL normdu  = 0.0;
    STATE k=1;
    for (int i = 1; i <= n_it; ++i) {
        // monta tangente e resíduo na solução atual x
        //std::cout << "Assembling.."   <<endl;
        Assemble();

        //std::cout << "Solving.."   <<endl;
        // resolve Δx
        Solve();
        dx = Solution();

        // atualiza solução candidata
        x += dx;
        LoadSolution(x);

        // *** RECOMPUTA o resíduo para a solução ATUALIZADA ***
        AssembleResidual();


        // métricas
        normdu  = Norm(dx);
        normrhs = Norm(Rhs());

        r0=r1;
        r1=r2;
        r2=normrhs;

        mr0=mr1;
        mr1=mr2;
        mr2=Rhs();


        if(verbose)
        {
            std::cout << " \n [it " << i << "] "
            << "||Δu|| = " << normdu
            << " | ||R|| = " << normrhs
            << " | tol = " << tol;// << std::endl;
            if (i > 3) {
                STATE lnR0 = log(r0), lnR1 = log(r1), lnR2 = log(r2);
                STATE p_est = (lnR2 - lnR1) / (lnR1 - lnR0);
                std::cout << " | p = " << p_est;
            }
        }



        if(normdu>10)return false;


        iters = i;

        // critério de convergência (resíduo relativo)
        if (normrhs < tol &&normdu<tol) {
            //std::cout << "normrhs ="<< normrhs<< " normdu ="<< normdu   <<endl;
            return true;
        }

    }

    std::cout << "NAO Convergiu com   ="<< iters<< " iteraçoes."<<std::endl;
    return false;


}
bool TPZElastoPlasticAnalysis::NewtonRaphson(REAL tol,TPZStack<STATE> &outresF,TPZStack<STATE> &outresU)
{

    TPZFMatrix<STATE> x(Solution()), dx(Solution());
    x.Zero(); dx.Zero();


    const int  n_it  = 30;
    const REAL EPS   = 1.e-30;

    //std::cout << "AssembleResidual.."   <<endl;
    int iters;
    AssembleResidual();
    REAL normrhs0 = Norm(fRhs);
    if (normrhs0 <1.e-3) normrhs0 = 1.; // proteção
    STATE r0=0.,r1=0.,r2=0.;
    TPZFMatrix<STATE> mr0,mr1,mr2;
    STATE rate0=0.,rate1=0.,rate2=0.;
    REAL normrhs = 1.0;
    REAL normdu  = 0.0;
    STATE k=1;
    for (int i = 1; i <= n_it; ++i) {
        // monta tangente e resíduo na solução atual x
        //std::cout << "Assembling.."   <<endl;
        Assemble();

        //std::cout << "Solving.."   <<endl;
        // resolve Δx
        Solve();
        dx = Solution();

        // atualiza solução candidata
        x += dx;
        LoadSolution(x);

        // *** RECOMPUTA o resíduo para a solução ATUALIZADA ***
        AssembleResidual();


        // métricas
        normdu  = Norm(dx);
        normrhs = Norm(Rhs());
        outresF.Push(normrhs);
        outresU.Push(normdu);
        r0=r1;
        r1=r2;
        r2=normrhs;

        mr0=mr1;
        mr1=mr2;
        mr2=Rhs();


        std::cout << " \n [it " << i << "] "
        << "||Δu|| = " << normdu
        << " | ||R|| = " << normrhs
        << " | tol = " << tol;// << std::endl;
        if (i > 3) {
            STATE lnR0 = log(r0), lnR1 = log(r1), lnR2 = log(r2);
            STATE p_est = (lnR2 - lnR1) / (lnR1 - lnR0);
            std::cout << " | p = " << p_est;
        }

        iters = i;

        // critério de convergência (resíduo relativo)
        if (normrhs < tol ) {
            std::cout << "normrhs ="<< normrhs<< " normdu ="<< normdu   <<endl;
            return true;
        }

    }

    std::cout << "NAO Convergiu com   ="<< iters<< " iteraçoes."<<std::endl;
    return false;


}

bool TPZElastoPlasticAnalysis::IterativeProcess(std::ostream &out,REAL tol, int numiter,bool linesearch, bool checkconv,int &iters)
{
    int    iter = 0;

    const int numeq = fCompMesh->NEquations();

    TPZFMatrix<STATE> prevsol(fSolution);
    if (prevsol.Rows() != numeq) prevsol.Redim(numeq,1);

    if (checkconv) {
        TPZVec<REAL> coefs(1,1.);
        TPZFMatrix<STATE> range(numeq,1,1.);
        CheckConvergence(*this, fSolution, range, coefs);
    }

    // resíduo inicial
    Assemble();
    REAL normrhs0 = Norm(fRhs);
    if (normrhs0 <1.e-3) normrhs0 = 1.; // proteção

    REAL normu = Norm(fSolution);
    REAL normu0 = Norm(fSolution);
    if (normu0  <1.e-3) normu0 = 1.; // proteção
    bool converged = false;
    // int maxiter=30;
    // if(numiter>maxiter) numiter=maxiter;

    // out << "\n[IterativeProcess2] Início do Newton-Raphson\n";
    // out << "  NumEq = " << numeq<< " | tol_u = " << tol << " | tol_f = " << tol<< " | maxiter = " << numiter << "\n";
     cout << "  Norma inicial do resíduo = " << normrhs0 << "\n";

    while (iter < numiter) {

        // monta e resolve o incremento
        Assemble();
        Solve(); // fSolution contém incremento (Δu) OU absol. (depende da sua infra)

        if (linesearch) {
            TPZFMatrix<STATE> nextSol;
            const REAL ls_tol = 0.1;
            const int  ls_it  = 200;

            switch (fLineSearch) {
                case ELineSearch::Armijo:
                    ArmijoLineSearch(prevsol, fSolution, nextSol, ls_tol, ls_it);
                    break;
                case ELineSearch::QuadraticArmijo:
                    QuadraticArmijoLineSearch(prevsol, fSolution, nextSol, ls_tol, ls_it);
                    break;
                case ELineSearch::GoldenSection:
                    GoldenSectionLineSearch(prevsol, fSolution, nextSol, ls_tol, ls_it);
                    break;
                case ELineSearch::Dicotomic:
                    DicotomicLineSearch(prevsol, fSolution, nextSol, ls_tol, ls_it);
                    break;
                case ELineSearch::StrongWolfe:
                    DebugStop();
                    StrongWolfeLineSearch(prevsol, fSolution, nextSol, ls_tol, ls_it);
                    break;
                case ELineSearch::NonmonotoneArmijo:
                    DebugStop();
                    NonmonotoneArmijoGLL(prevsol,fSolution,nextSol,fPhiHistory,1e-4,0.5,7,1.0,ls_it);
                    break;
                case ELineSearch::None:
                default:
                    break;
            }
            fSolution = nextSol;
            //out << "  [it " << iter << "] line search aplicado\n";
        } else {
            fSolution += prevsol;
        }

        // erro de deslocamento
        TPZFMatrix<STATE> delta = fSolution;
        delta -= prevsol;
        const REAL err_u = Norm(delta);

        // atualiza estado interno
        prevsol = fSolution;
        LoadSolution(fSolution);

        // reavalia resíduo no estado ATUAL
        AssembleResidual();
        const REAL err_f = Norm(fRhs) ;

        // imprime diagnóstico
        std::cout << "  [it " << iter << "] "
        << "||Δu|| = " << err_u
        << " | ||R|| = " << err_f
        << " | tol = " << tol << std::endl;

        // critério de parada
        if ( err_u<tol && err_f<tol*100)
        //if ( err_u<tol)
        {
            out << "  -> Convergência atingida em " << iter+1 << " iterações " << " ||Δu||/||Δu0|| = " << err_u<< " | ||R||/||R0|| = " << err_f<< " | tol = " << tol << "\n";
            converged = true;
            iter++;
            break;
        }

        iter++;
    }

    iters = iter;
    if(!converged){
        //out << "  -> Não convergiu em " << iters<< " iterações. Últimos erros: "<< "||Δu||/||Δu0||=" << prev_err_u<< " | ||R||/||R0||=" << prev_err_f << "\n";
    }
    out.flush();
    return converged;
}

void TPZElastoPlasticAnalysis::TransferSolution()
{

}
REAL TPZElastoPlasticAnalysis::NonmonotoneArmijoGLL(const TPZFMatrix<STATE>& Wn,
                          const TPZFMatrix<STATE>& d,
                          TPZFMatrix<STATE>& NextW,
                          std::deque<REAL>& phi_hist, // mantém últimas M φ
                          REAL c1, REAL beta,
                          int M, REAL a0, int max_red)
{
    // φ(0) e φ'(0)
    TPZFMatrix<REAL> Rn; this->Residual(Rn, 0);
    REAL phi0 = (REAL)0.5 * Dot(Rn, Rn);
    REAL dphi0 = -Dot(Rn, Rn); // se Newton; caso contrário compute Jd como no Wolfe
    REAL PhiMax = phi0;
    for (REAL v : phi_hist) PhiMax = std::max(PhiMax, v);

    REAL a = a0;
    for (int k=0; k<max_red; ++k) {
        TPZFMatrix<STATE> Wtry = Wn; Wtry += a*d;
        TPZFMatrix<REAL> Rtry; this->Residual(Rtry, 0);
        REAL phia = (REAL)0.5 * Dot(Rtry, Rtry);
        if (phia <= PhiMax + c1*a*dphi0) { // aceita
            NextW = Wtry;
            phi_hist.push_back(phia);
            if ((int)phi_hist.size() > M) phi_hist.pop_front();
            return a;
        }
        // backtracking (pode trocar por interpolação quadrática protegida)
        a *= beta;
    }
    NextW = Wn; NextW += a*d;
    return a;
}

REAL TPZElastoPlasticAnalysis::StrongWolfeLineSearch(const TPZFMatrix<STATE>& Wn,const TPZFMatrix<STATE>& d,TPZFMatrix<STATE>& NextW,REAL c1, REAL c2,REAL a0, int max_eval)
{
    auto phi = [&](const TPZFMatrix<STATE>& W, TPZFMatrix<STATE>& R)->REAL{
        TPZFMatrix<STATE> Wtmp(W); this->LoadSolution(Wtmp);
        TPZFMatrix<REAL> res; this->Residual(res, 0);
        return (REAL)0.5 * Dot(res, res); // Inner = R^T R
    };
    auto dphi = [&](const TPZFMatrix<STATE>& W,
                    const TPZFMatrix<STATE>& R)->REAL{
        // monta J(W)*d -> Jd
        TPZFMatrix<REAL> K; TPZVec<REAL> coefs(1,1.0); this->ComputeTangent(K, coefs, 0);
        TPZFMatrix<REAL> Jd; K.Multiply(d, Jd);
        return Dot(R, Jd); // R^T (J d)
    };

    TPZFMatrix<STATE> Rn; this->Residual(Rn, 0);
    REAL phi0 = (REAL)0.5 * Dot(Rn, Rn);
    REAL dphi0;
    { // se for Newton puro, use atalho; caso contrário, compute
      // dphi0 = -||R||^2;
      TPZFMatrix<REAL> K; TPZVec<REAL> coefs(1,1.0); this->ComputeTangent(K, coefs, 0);
      TPZFMatrix<REAL> Jd; K.Multiply(d, Jd);
      dphi0 = Dot(Rn, Jd);
    }

    REAL alo=0, ahi=a0, philo=phi0, dphilo=dphi0;
    TPZFMatrix<STATE> Wtrial, Rtrial;
    for (int k=0; k<max_eval; ++k) {
        // avalia em ahi
        Wtrial = Wn; Wtrial += ahi * d;
        REAL phihi = phi(Wtrial, Rtrial);
        if ( (phihi > phi0 + c1*ahi*dphi0) || (k>0 && phihi >= philo) ) {
            // entra no "zoom"
            REAL aL=alo, aH=ahi; REAL phiL=philo, dphiL=dphilo;
            for (int z=0; z<max_eval; ++z) {
                // interpolação cúbica protegida entre [aL,aH]
                REAL a = 0.5*(aL+aH);
                TPZFMatrix<STATE> Wz = Wn; Wz += a*d;
                TPZFMatrix<STATE> Rz; REAL phiz = phi(Wz, Rz);
                if ( (phiz > phi0 + c1*a*dphi0) || (phiz >= phiL) ) {
                    aH = a;
                } else {
                    REAL dphiz = dphi(Wz, Rz);
                    if ( std::fabs(dphiz) <= c2*std::fabs(dphi0) ) {
                        NextW = Wz; return a;
                    }
                    if ( (aH - aL)*dphiz >= 0 ) aH = aL;
                    aL = a; phiL = phiz; dphiL = dphiz;
                }
                if (std::fabs(aH-aL) < 1e-12) { NextW = Wz; return a; }
            }
        }
        REAL dphihi = dphi(Wtrial, Rtrial);
        if ( std::fabs(dphihi) <= c2*std::fabs(dphi0) ) { NextW = Wtrial; return ahi; }
        if ( dphihi >= 0 ) {
            // entra no "zoom" com bracket [ahi, alo]
            REAL aL=ahi, aH=alo; std::swap(aL,aH); // garanta aL<->alo
            // (mesma rotina de zoom acima…)
        }
        // expande
        alo = ahi; philo = phihi; dphilo = dphihi; ahi *= 2.0;
    }
    // falha branda: devolve o melhor visto
    NextW = Wn; NextW += alo * d; return alo;
}

REAL TPZElastoPlasticAnalysis::DicotomicLineSearch(const TPZFMatrix<STATE>& Wn,
                                                   TPZFMatrix<STATE> DeltaW,
                                                   TPZFMatrix<STATE>& NextW,
                                                   REAL tol, int niter)
{
    // intervalo inicial de α
    REAL A = (REAL)0.0, B = (REAL)1.0;

    const REAL delta_min = (REAL)1e-6;
    auto delta_for = [&](REAL width){
        return std::max<REAL>(delta_min, (REAL)0.1 * width); // 10% da largura
    };

    const REAL amin = std::max<REAL>(tol, (REAL)1e-8); // alpha mínimo aceitável

    // valida DeltaW
    if (DeltaW.Rows() == 0 || DeltaW.Cols() == 0) {
        NextW = Wn;
        // aplica um pequeno passo para evitar alpha == 0 problema downstream
        TPZFMatrix<STATE> tiny = Wn;
        TPZFMatrix<STATE> dd = DeltaW;
        dd *= amin;
        tiny += dd;
        this->LoadSolution(tiny);
        NextW = tiny;
        return amin;
    }

    // --- preparar análise temporária CLONANDO a malha e materiais ---
    TPZCompMesh *origMesh = this->Mesh();
    if (!origMesh) {
        throw std::runtime_error("DicotomicLineSearch: mesh nula no objeto this.");
    }

    TPZCompMesh *meshCopy = nullptr;
    try {
        // Substitua Clone() pelo método correto da sua versão do NeoPZ se necessário.
        meshCopy = origMesh->Clone();
    } catch (...) {
        meshCopy = nullptr;
    }

    if (!meshCopy) {
        throw std::runtime_error("DicotomicLineSearch: clonagem profunda da malha nao disponivel. "
        "Implemente clonagem de estados internos ou use snapshot.");
    }

    // Constrói uma análise temporária com a mesh copiada.
    // Assumimos que tmpAnalysis gerencia a mesh copiada (ajuste se sua API for diferente).
    TPZElastoPlasticAnalysis tmpAnalysis(meshCopy,std::cout);
    // Opcional: copie solver / structural matrix / configurações relevantes do this para tmpAnalysis
    // tmpAnalysis.SetStructuralMatrix(this->StructuralMatrix()); // adapte conforme API
    // tmpAnalysis.SetSolver(this->Solver()); // adapte conforme API

    // Função objetivo sem efeitos colaterais sobre `this` (avalia em tmpAnalysis)
    auto eval_phi_tmp = [&](const TPZFMatrix<STATE>& Wc)->REAL {
        TPZFMatrix<STATE> W = Wc;           // cópia para LoadSolution
        tmpAnalysis.LoadSolution(W);
        tmpAnalysis.AssembleResidual();
        REAL nR = Norm(tmpAnalysis.fRhs);
        if (!std::isfinite(nR)) return std::numeric_limits<REAL>::infinity();
        // Usar objetivo 0.5 * ||R||^2 para coerência com line-search clássico
        return (REAL)0.5 * nR * nR;
    };

    // computa phi0 em Wn (usando tmpAnalysis para não alterar this)
    const REAL phi0 = eval_phi_tmp(Wn);

    // iteração dicotômica
    int it = 0;
    REAL width = B - A;

    REAL best_phi = std::numeric_limits<REAL>::infinity();
    REAL best_alpha = (REAL)0.5 * (A + B);
    TPZFMatrix<STATE> bestW = Wn;

    while (it < niter && width > tol) {
        const REAL mid = (REAL)0.5*(A + B);
        REAL delta = delta_for(width);

        REAL x1 = std::max<REAL>(A, mid - delta);
        REAL x2 = std::min<REAL>(B, mid + delta);
        if (x1 >= x2) {
            // relaxa delta simetricamente; se ainda não houver espaço, sai
            delta *= (REAL)0.5;
            x1 = std::max<REAL>(A, mid - delta);
            x2 = std::min<REAL>(B, mid + delta);
            if (x1 >= x2) break;
        }

        TPZFMatrix<STATE> t1 = Wn, t2 = Wn;
        TPZFMatrix<STATE> d1 = DeltaW, d2 = DeltaW;
        d1 *= x1; t1 += d1;
        d2 *= x2; t2 += d2;

        const REAL f1 = eval_phi_tmp(t1);
        const REAL f2 = eval_phi_tmp(t2);

        if (std::isnan(f1) || std::isnan(f2)) {
            // avaliação inválida: interrompe busca (ponto de emergência)
            break;
        }

        if (f1 < best_phi) { best_phi = f1; best_alpha = x1; bestW = t1; }
        if (f2 < best_phi) { best_phi = f2; best_alpha = x2; bestW = t2; }

        // dicotômico: manter metade com menor valor
        if (f1 > f2) {
            A = x1; // mínimo em (x1,B]
        } else if (f2 > f1) {
            B = x2; // mínimo em [A,x2)
        } else {
            // empate: reduzir simetricamente para evitar viés
            A = x1;
            B = x2;
        }

        width = B - A;
        ++it;
    }

    // escolha final: prefira melhor amostrado; se não encontrado, meio do intervalo
    REAL alpha = (best_phi < std::numeric_limits<REAL>::infinity()) ? best_alpha : (REAL)0.5*(A + B);
    if (alpha < amin) alpha = amin;

    // monta NextW e COMITA na análise real (aplica o passo)
    NextW = Wn;
    TPZFMatrix<STATE> DeltaScaled = DeltaW; // cópia segura
    DeltaScaled *= alpha;
    NextW += DeltaScaled;

    // Aplica/commit NextW na análise real para que o próximo solver trabalhe com o novo estado
    this->LoadSolution(NextW);

    // NOTA: não deletamos meshCopy explicitamente; espera-se que tmpAnalysis libere a mesh copiada em seu destrutor.
    // Se sua API exigir delete(meshCopy), adapte aqui.

    return alpha;
}


REAL TPZElastoPlasticAnalysis::QuadraticArmijoLineSearch(const TPZFMatrix<STATE>& Wn,
                                                         TPZFMatrix<STATE> DeltaW,
                                                         TPZFMatrix<STATE>& NextW,
                                                         REAL tol, int niter)
{
    const REAL c1        = (REAL)1e-4;                // Armijo parameter
    const REAL amin      = std::max<REAL>(tol, (REAL)1e-8);
    const REAL shrink_lo = (REAL)0.1, shrink_hi = (REAL)0.5; // proteção da interpolação

    // utilitário: φ(W) = 0.5 ||R||^2 (avalia em cópia para LoadSolution)
    auto eval_phi = [&](const TPZFMatrix<STATE>& Wc)->REAL {
        TPZFMatrix<STATE> W = Wc;               // NÃO-const p/ LoadSolution
        this->LoadSolution(W);
        this->AssembleResidual();               // só o resíduo
        const REAL rn = Norm(this->fRhs);
        if (!std::isfinite(rn)) return std::numeric_limits<REAL>::infinity();
        return (REAL)0.5 * rn * rn;
    };

    // Guard RAII para restaurar fSolution ao sair, a menos que Commit seja chamado.
    TPZFMatrix<STATE> backup = this->fSolution;
    struct RestoreGuard {
        TPZElastoPlasticAnalysis* an;
        TPZFMatrix<STATE> backup;
        bool committed;
        RestoreGuard(TPZElastoPlasticAnalysis* a, const TPZFMatrix<STATE>& b)
        : an(a), backup(b), committed(false) {}
        ~RestoreGuard() {
            if (!committed && an) {
                an->LoadSolution(backup);
            }
        }
        void Commit() { committed = true; }
    } guard(this, backup);

    // φ(0) e φ'(0) aproximado
    TPZFMatrix<STATE> W0 = Wn;                  // cópia não-const
    this->LoadSolution(W0);
    this->AssembleResidual();
    const REAL r0 = Norm(this->fRhs);
    REAL phi0 = (REAL)0.5 * r0 * r0;
    // aproximação inicial para φ'(0). No contexto do resíduo, usa-se -||R||^2 como heurística.
    REAL phip0 = -(REAL)(r0 * r0);

    // fallback para phip0 se não for descida
    if (phip0 >= (REAL)0) {
        const REAL eps = (REAL)1e-6;
        TPZFMatrix<STATE> Wp = Wn;
        TPZFMatrix<STATE> dd = DeltaW;
        dd *= eps;
        Wp += dd;
        const REAL phip = eval_phi(Wp);
        if (std::isfinite(phip)) {
            phip0 = (phip - phi0) / eps;
        } else {
            phip0 = -(REAL)std::max((REAL)1e-16, phi0);
        }
        if (phip0 >= (REAL)0) phip0 = -(REAL)std::max((REAL)1e-16, phi0);
    }

    // se ΔW vazio: não anda (restauração automática pelo guard). retorna amin para evitar 0.
    if (DeltaW.Rows()==0 || DeltaW.Cols()==0) {
        NextW = Wn;
        return amin;
    }

    REAL alpha = (REAL)1.0;
    REAL best_phi = phi0, best_alpha = amin;
    TPZFMatrix<STATE> bestW = Wn;

    for (int k = 0; k < niter && alpha >= amin; ++k) {
        TPZFMatrix<STATE> trial = Wn;
        TPZFMatrix<STATE> d = DeltaW;
        d *= alpha;
        trial += d;

        const REAL phi_a = eval_phi(trial);
        if (!std::isfinite(phi_a)) {
            // avaliação inválida: reduzir e continuar (proteção)
            alpha *= shrink_hi;
            continue;
        }

        if (phi_a < best_phi) { best_phi = phi_a; best_alpha = alpha; bestW = trial; }

        // Armijo: φ(α) ≤ φ(0) + c1 * α * φ'(0)
        if (phi_a <= (REAL)(phi0 + c1 * alpha * phip0)) {
            // aceita: aplica trial e comita (não restaurar)
            this->LoadSolution(trial);
            guard.Commit();
            NextW = trial;
            return std::max<REAL>(alpha, amin);
        }

        // --- interpolação quadrática protegida ---
        // modelo: φ(α) ≈ φ0 + φ'0 α + c α^2 => c = (φ(α) - φ0 - φ'0 α)/α^2
        REAL denom = (phi_a - phi0 - phip0 * alpha);
        REAL a_quad;
        if (denom <= (REAL)0) {
            // denom não positivo => não confiamos na interpolação, reduzimos "safely"
            a_quad = alpha * shrink_hi;
        } else {
            a_quad = -(phip0) * alpha * alpha / ((REAL)2.0 * denom);
            // proteger dentro de [shrink_lo*alpha, shrink_hi*alpha]
            a_quad = std::max(shrink_lo * alpha, std::min(shrink_hi * alpha, a_quad));
        }
        // garantir que alpha diminua pra evitar loop infinito (proteção extra)
        if (a_quad >= alpha) {
            alpha *= shrink_hi;
        } else {
            alpha = a_quad;
        }
    }

    // fallback: se alguma amostra melhorou, aplique o melhor; senão aplique um tiny step amin
    if (best_phi < phi0) {
        this->LoadSolution(bestW);
        guard.Commit();
        NextW = bestW;
        return std::max<REAL>(best_alpha, amin);
    } else {
        // aplica tiny step para evitar Δu==0 no passo seguinte
        TPZFMatrix<STATE> tiny = Wn;
        TPZFMatrix<STATE> dd = DeltaW;
        dd *= amin;
        tiny += dd;
        this->LoadSolution(tiny);
        guard.Commit();
        NextW = tiny;
        return amin;
    }
}
// Armijo (backtracking) line search — evita retornar alpha == 0 e garante que,
// se um trial for aceito, a solução é aplicada (committed) na análise.
// Avaliações são feitas no próprio `this` com backup/restore; AO ACEITAR,
// o guard é marcado como committed para não restaurar o estado antigo.
//
// Observações:
// - Se a sua versão do NeoPZ suporta clonagem profunda (TPZCompMesh::Clone() etc.),
//   é preferível avaliar em uma análise clonada e no fim aplicar LoadSolution(NextW)
//   em `this`. Aqui fazemos a versão que compila com LoadSolution(TPZFMatrix<STATE>&).

REAL TPZElastoPlasticAnalysis::ArmijoLineSearch(const TPZFMatrix<STATE>& Wn,
                                                TPZFMatrix<STATE> DeltaW,
                                                TPZFMatrix<STATE>& NextW,
                                                REAL tol, int niter)
{
    const REAL c    = (REAL)1e-4;                   // parâmetro Armijo
    const REAL rho  = (REAL)0.5;                    // redução do passo
    const REAL amin = std::max<REAL>(tol, (REAL)1e-8); // alpha mínimo aceitável

    // Guard RAII para restaurar fSolution ao sair, a menos que commit() seja chamado.
    TPZFMatrix<STATE> backup = this->fSolution;
    struct RestoreGuard {
        TPZElastoPlasticAnalysis* an;
        TPZFMatrix<STATE> backup;
        bool committed;
        RestoreGuard(TPZElastoPlasticAnalysis* a, const TPZFMatrix<STATE>& b)
        : an(a), backup(b), committed(false) {}
        ~RestoreGuard() {
            if (!committed && an) {
                an->LoadSolution(backup);
            }
        }
        void Commit() { committed = true; }
    } guard(this, backup);

    // objetivo: f = 0.5 * ||R||^2 (mais coerente com Armijo)
    auto eval_f = [&](const TPZFMatrix<STATE>& W)->REAL {
        TPZFMatrix<STATE> tmp = W;        // cópia NÃO-CONST para LoadSolution
        this->LoadSolution(tmp);
        this->AssembleResidual();         // monta apenas o resíduo
        REAL nR = Norm(this->fRhs);
        if (!std::isfinite(nR)) return std::numeric_limits<REAL>::infinity();
        return (REAL)0.5 * nR * nR;
    };

    // valida DeltaW
    if (DeltaW.Rows() == 0 || DeltaW.Cols() == 0) {
        NextW = Wn;
        // não commit: restauração automática no guard
        return amin;
    }

    // valor inicial f(0)
    const REAL phi0 = eval_f(Wn);

    // inicial
    REAL alpha = (REAL)1.0;

    // melhor já visto (fallback)
    REAL best_phi = phi0;
    REAL best_alpha = amin; // inicializamos com amin para evitar 0
    TPZFMatrix<STATE> bestW = Wn;

    for (int k = 0; k < niter && alpha >= amin; ++k) {
        TPZFMatrix<STATE> trial = Wn;
        TPZFMatrix<STATE> d = DeltaW;
        d *= alpha;
        trial += d;

        const REAL phi_a = eval_f(trial);
        if (!std::isfinite(phi_a)) {
            // avaliação inválida: reduzir e continuar
            alpha *= rho;
            continue;
        }

        if (phi_a < best_phi) { best_phi = phi_a; best_alpha = alpha; bestW = trial; }

        // critério Armijo adaptado (usando f): f(alpha) <= (1 - c*alpha) * f0
        if (phi_a <= (REAL)((1.0 - c * alpha) * phi0)) {
            // aceita: aplica trial na análise real e comita (não restaurar)
            this->LoadSolution(trial);   // aplicar a solução aceita
            guard.Commit();              // impede restauração no destrutor
            NextW = trial;
            return std::max<REAL>(alpha, amin);
        }

        alpha *= rho;
    }

    // fallback: se alguma amostra melhorou, aplique-a; caso contrário aplique um passo mínimo
    if (best_phi < phi0) {
        this->LoadSolution(bestW);
        guard.Commit();
        NextW = bestW;
        return std::max<REAL>(best_alpha, amin);
    } else {
        // não houve melhoria: aplique um pequeno passo amin (evita Δu == 0)
        TPZFMatrix<STATE> tiny = Wn;
        TPZFMatrix<STATE> dd = DeltaW;
        dd *= amin;
        tiny += dd;
        this->LoadSolution(tiny);
        guard.Commit();
        NextW = tiny;
        return amin;
    }
}

// Golden-section line search — versão corrigida e robusta.
// - Avaliações de φ feitas em uma análise temporária clonada (tmpAnalysis) para NÃO alterar `this`.
// - Objetivo usado: f = 0.5 * ||R||^2 (coerente com armijo/quadratic methods).
// - Nunca retorna alpha == 0: impõe amin = max(tol, 1e-8) e aplica tiny step se necessário.
// - Aplica (commit) NextW na análise real (this->LoadSolution) ao final.
// - Trata avaliações inválidas (NaN/Inf) e escolhe o melhor ponto amostrado como fallback.
//
// Observações:
// - Presumo existência de TPZCompMesh::Clone() e um construtor TPZElastoPlasticAnalysis(TPZCompMesh*).
//   Se a sua API for diferente, adapte as chamadas de clonagem/construct conforme necessário.
// - Aqui assumo que tmpAnalysis assume a propriedade da mesh copiada; não faço delete(meshCopy).
//   Se sua API exigir explicitamente liberar meshCopy, ajuste o código.

REAL TPZElastoPlasticAnalysis::GoldenSectionLineSearch(const TPZFMatrix<STATE>& Wn,
                                                       TPZFMatrix<STATE> DeltaW,
                                                       TPZFMatrix<STATE>& NextW,
                                                       REAL tol, int niter)
{
    // parámetros
    const REAL amin = std::max<REAL>(tol, (REAL)1e-8); // alpha mínimo aceitável
    constexpr bool kVerbose = false;

    // valida DeltaW
    if (DeltaW.Rows() == 0 || DeltaW.Cols() == 0) {
        // aplica tiny step e comita para evitar alpha==0 downstream
        NextW = Wn;
        TPZFMatrix<STATE> tiny = DeltaW;
        tiny *= amin;
        NextW += tiny;
        this->LoadSolution(NextW);
        return amin;
    }

    // --- preparar análise temporária CLONANDO a malha e materiais ---
    TPZCompMesh *origMesh = this->Mesh();
    if (!origMesh) {
        throw std::runtime_error("GoldenSectionLineSearch: mesh nula no objeto this.");
    }

    TPZCompMesh *meshCopy = nullptr;
    try {
        meshCopy = origMesh->Clone();
    } catch (...) {
        meshCopy = nullptr;
    }

    if (!meshCopy) {
        throw std::runtime_error("GoldenSectionLineSearch: clonagem profunda da malha nao disponivel. "
        "Implemente clonagem de estados internos ou use snapshot.");
    }

    // Constrói uma análise temporária com a mesh copiada.
    TPZElastoPlasticAnalysis tmpAnalysis(meshCopy,std::cout);
    // opcional: copiar configurações do solver/structural matrix se necessário
    // tmpAnalysis.SetStructuralMatrix(this->StructuralMatrix()); // adaptar conforme API
    // tmpAnalysis.SetSolver(this->Solver()); // adaptar conforme API

    // utilitário: avalia φ = 0.5 * ||R||^2 em tmpAnalysis (não altera `this`)
    auto EvalPhiTmp = [&](REAL alpha)->REAL {
        TPZFMatrix<STATE> trial = Wn;
        TPZFMatrix<STATE> d = DeltaW;
        d *= alpha;
        trial += d;
        tmpAnalysis.LoadSolution(trial);
        tmpAnalysis.AssembleResidual();
        REAL nr = Norm(tmpAnalysis.fRhs);
        if (!std::isfinite(nr)) return std::numeric_limits<REAL>::infinity();
        return (REAL)0.5 * nr * nr;
    };

    // Golden-section constants
    constexpr REAL gr = (REAL)0.6180339887498949;  // phi
    constexpr REAL gr2 = (REAL)1.0 - gr;           // 0.381966...

    // extremos
    REAL A = (REAL)0.0;
    REAL B = (REAL)1.0;

    REAL f0 = EvalPhiTmp((REAL)0.0);
    REAL f1 = EvalPhiTmp((REAL)1.0);

    // interior points
    REAL L = A + gr2 * (B - A);
    REAL M = A + gr  * (B - A);

    REAL fL = EvalPhiTmp(L);
    REAL fM = EvalPhiTmp(M);

    int it = 0;
    REAL width = B - A;
    int res_evals = 4;

    // track best sampled
    REAL best_phi = std::numeric_limits<REAL>::infinity();
    REAL best_alpha = (REAL)0.5 * (A + B);
    TPZFMatrix<STATE> bestW = Wn;

    auto consider_sample = [&](REAL alpha, REAL phi){
        if (std::isfinite(phi) && phi < best_phi) {
            best_phi = phi;
            best_alpha = alpha;
            // build bestW lazily
            bestW = Wn;
            TPZFMatrix<STATE> dd = DeltaW;
            dd *= alpha;
            bestW += dd;
        }
    };

    consider_sample((REAL)0.0, f0);
    consider_sample((REAL)1.0, f1);
    consider_sample(L, fL);
    consider_sample(M, fM);

    if (kVerbose) {
        std::cout << "[Golden] start f0="<<f0<<" f1="<<f1<<" L="<<L<<" fL="<<fL<<" M="<<M<<" fM="<<fM<<"\n";
    }

    // loop golden
    while (it < niter && width > tol) {
        if (fL > fM) {
            // minimum in (L, B]
            A = L;
            L = M;
            fL = fM;
            M = A + gr * (B - A);
            fM = EvalPhiTmp(M);
            ++res_evals;
            consider_sample(M, fM);
        } else {
            // minimum in [A, M)
            B = M;
            M = L;
            fM = fL;
            L = A + gr2 * (B - A);
            fL = EvalPhiTmp(L);
            ++res_evals;
            consider_sample(L, fL);
        }
        width = B - A;
        ++it;
        if (kVerbose) {
            std::cout << " [it " << it << "] A="<<A<<" B="<<B<<" L="<<L<<" fL="<<fL<<" M="<<M<<" fM="<<fM<<" width="<<width<<"\n";
        }
        // protection: if both fL and fM are infinite/NaN, abort
        if (!std::isfinite(fL) && !std::isfinite(fM)) break;
    }

    // choose alpha: prefer best sampled; otherwise midpoint
    REAL alpha = (best_phi < std::numeric_limits<REAL>::infinity()) ? best_alpha : (REAL)0.5*(A + B);
    if (alpha < amin) alpha = amin;

    // build NextW and commit to real analysis
    NextW = Wn;
    TPZFMatrix<STATE> dsc = DeltaW;
    dsc *= alpha;
    NextW += dsc;

    // Apply the chosen solution to the real analysis so the next solver sees it
    this->LoadSolution(NextW);

    if (kVerbose) {
        std::cout << "[Golden] finish alpha=" << alpha
        << " iters=" << it
        << " evals=" << res_evals
        << " best_phi=" << best_phi << "\n";
    }

    return alpha;
}
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
        cout<< "the material is linear, exiting..."<<endl;
        return 0.;
    }


	if(ResetOutputDisplacements)
	{
		fCumSol.Zero();
	}else{
        //cout<< "accumulating solution..."<<endl;
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



void TPZElastoPlasticAnalysis::SetAllCreateFunctionsWithMem(TPZCompMesh *cmesh)
{
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
