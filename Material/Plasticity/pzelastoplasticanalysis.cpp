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
        dx = Solution();
        x += dx;
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


bool TPZElastoPlasticAnalysis::IterativeProcess(std::ostream &out,REAL tol, int numiter,bool linesearch, bool checkconv,int &iters)
{
    int    iter = 0;
    REAL   prev_err_u = std::numeric_limits<REAL>::infinity();
    REAL   prev_err_f = std::numeric_limits<REAL>::infinity();



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
    if (normrhs0 == 0.) normrhs0 = 1.; // proteção

    REAL normu0 = Norm(fSolution);
    if (normu0 == 0.) normu0 = 1.; // proteção
    bool converged = false;
    // int maxiter=30;
    // if(numiter>maxiter) numiter=maxiter;

    // out << "\n[IterativeProcess2] Início do Newton-Raphson\n";
    // out << "  NumEq = " << numeq<< " | tol_u = " << tol << " | tol_f = " << tol<< " | maxiter = " << numiter << "\n";
    // out << "  Norma inicial do resíduo = " << normrhs0 << "\n";

    while (iter < numiter) {

        // monta e resolve o incremento
        Assemble();
        Solve(); // fSolution contém incremento (Δu) OU absol. (depende da sua infra)

        if (linesearch) {
            TPZFMatrix<STATE> nextSol;
            const REAL ls_tol = (REAL)1e-3 * std::max<REAL>( (REAL)1.0, Norm(fSolution) );
            const int  ls_it  = 60;
            if(true)
            {
                ArmijoLineSearch(prevsol, fSolution, nextSol, ls_tol, ls_it); // nextSol = solução ABSOLUTA
                //QuadraticArmijoLineSearch(prevsol, fSolution, nextSol, ls_tol, ls_it);
            }else{
                DicotomicLineSearch(prevsol, fSolution, nextSol, ls_tol, ls_it); // nextSol = solução ABSOLUTA
            }

            fSolution = nextSol;
            //out << "  [it " << iter << "] line search aplicado\n";
        } else {
            fSolution += prevsol;
        }

        // erro de deslocamento
        TPZFMatrix<STATE> delta = fSolution;
        delta -= prevsol;
        const REAL err_u = Norm(delta)/normu0;

        // atualiza estado interno
        prevsol = fSolution;
        LoadSolution(fSolution);

        // reavalia resíduo no estado ATUAL
        AssembleResidual();
        const REAL err_f = Norm(fRhs) / normrhs0;

        // imprime diagnóstico
        // out << "  [it " << iter << "] "
        // << "||Δu||/normu0 = " << err_u/normu0
        // << " | ||R||/||R0|| = " << err_f
        // << " | tol = " << tol << "\n";

        // critério de parada
        //if (err_u <= tol && err_f <= tol) {
            if ( err_u <= tol && err_f <= tol) {
            //out << "  -> Convergência atingida em " << iter+1 << " iterações " << " ||Δu||/||Δu0|| = " << err_u<< " | ||R||/||R0|| = " << err_f<< " | tol = " << tol << "\n";
            converged = true;
            iter++;
            break;
        }

        prev_err_u = err_u;
        prev_err_f = err_f;

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
// Busca linear dicotômica (dichotomous search) para minimizar ||R|| ao longo de ΔW.
// tol   -> tolerância para a largura do intervalo em α (ex.: 1e-3)
// niter -> número máx. de iterações
REAL TPZElastoPlasticAnalysis::DicotomicLineSearch(const TPZFMatrix<STATE>& Wn,
                                          TPZFMatrix<STATE> DeltaW,
                                          TPZFMatrix<STATE>& NextW,
                                          REAL tol, int niter)
{
    // intervalo inicial de α
    REAL A = (REAL)0.0, B = (REAL)1.0;

    // pequeno deslocamento interno para as duas amostras (padrão da dicotomização)
    // mantém proporcionalidade ao intervalo, com piso para não colapsar
    const REAL delta_min = (REAL)1e-6;
    auto delta_for = [&](REAL width){
        return std::max<REAL>(delta_min, (REAL)0.1 * width); // 10% da largura
    };

    // utilitário: avalia ||R|| em um W qualquer (cópia não-const para LoadSolution)
    auto eval_phi = [&](const TPZFMatrix<STATE>& Wc)->REAL {
        TPZFMatrix<STATE> W = Wc;   // cópia NÃO-CONST
        this->LoadSolution(W);
        this->AssembleResidual();   // monta APENAS o resíduo
        return Norm(this->fRhs);    // ||R||_2
    };

    // backup do estado para restaurar no retorno
    TPZFMatrix<STATE> backup = this->fSolution;

    // se ΔW vazio, retorna sem andar
    if (DeltaW.Rows()==0 || DeltaW.Cols()==0) {
        NextW = Wn;
        this->LoadSolution(backup);
        return (REAL)0.0;
    }

    int  it    = 0;
    REAL width = B - A;

    // opcional: guarda melhor ponto já visto (robustez)
    REAL best_phi   = std::numeric_limits<REAL>::infinity();
    REAL best_alpha = (REAL)0.0;
    TPZFMatrix<STATE> bestW = Wn;

    while (it < niter && width > tol) {
        const REAL mid = (REAL)0.5*(A + B);
        REAL delta = delta_for(width);

        // garante que x1 < x2 e ambos dentro de [A,B]
        REAL x1 = std::max<REAL>(A, mid - delta);
        REAL x2 = std::min<REAL>(B, mid + delta);
        if (x1 >= x2) { // se encostou, relaxa delta
            delta *= (REAL)0.5;
            x1 = std::max<REAL>(A, mid - delta);
            x2 = std::min<REAL>(B, mid + delta);
            if (x1 >= x2) break; // não há espaço útil
        }

        // avalia φ nos dois pontos
        TPZFMatrix<STATE> trial1 = Wn, trial2 = Wn;
        TPZFMatrix<STATE> d1 = DeltaW, d2 = DeltaW;
        d1 *= x1; trial1 += d1;
        d2 *= x2; trial2 += d2;

        const REAL f1 = eval_phi(trial1);
        const REAL f2 = eval_phi(trial2);

        if (f1 < best_phi) { best_phi = f1; best_alpha = x1; bestW = trial1; }
        if (f2 < best_phi) { best_phi = f2; best_alpha = x2; bestW = trial2; }

        // passo dicotômico: mantém a metade com menor valor
        if (f1 > f2) {
            // mínimo está em (x1, B]
            A = x1;
        } else {
            // mínimo está em [A, x2)
            B = x2;
        }

        width = B - A;
        it++;
    }

    // α escolhido: meio do intervalo final; se quiser, pode usar best_alpha
    const REAL alpha = (REAL)0.5*(A + B);

    // monta NextW = Wn + α ΔW
    NextW = Wn;
    DeltaW *= alpha;  // DeltaW é por valor → seguro
    NextW += DeltaW;

    // restaura estado interno do Analysis
    this->LoadSolution(backup);

    return alpha;
}

// Armijo com interpolação quadrática protegida (rápido/robusto)
REAL TPZElastoPlasticAnalysis::QuadraticArmijoLineSearch(const TPZFMatrix<STATE>& Wn,
                                          TPZFMatrix<STATE> DeltaW,
                                          TPZFMatrix<STATE>& NextW,
                                          REAL tol, int niter)
{
    const REAL c1   = (REAL)1e-4;               // decréscimo suficiente
    const REAL amin = std::max<REAL>(tol, (REAL)1e-8);
    const REAL shrink_lo = (REAL)0.1, shrink_hi = (REAL)0.5; // guarda-corpos

    // utilitário: φ(W) = 1/2 ||R||^2
    auto eval_phi = [&](const TPZFMatrix<STATE>& Wc)->REAL {
        TPZFMatrix<STATE> W = Wc;               // NÃO-const p/ LoadSolution
        this->LoadSolution(W);
        this->AssembleResidual();               // só o resíduo
        const REAL rn = Norm(this->fRhs);
        return (REAL)0.5 * rn * rn;
    };

    // backup do estado interno
    TPZFMatrix<STATE> backup = this->fSolution;

    // φ(0) e φ'(0)
    TPZFMatrix<STATE> W0 = Wn;                  // cópia não-const
    this->LoadSolution(W0);
    this->AssembleResidual();
    const REAL r0  = Norm(this->fRhs);
    REAL phi0      = (REAL)0.5 * r0 * r0;
    REAL phip0     = -(REAL)(r0 * r0);          // para passo de Newton, φ'(0)= -||R||^2

    // fallback se, por algum motivo, φ'(0)≥0
    if (phip0 >= (REAL)0) {
        const REAL eps = (REAL)1e-6;
        TPZFMatrix<STATE> Wp = Wn, d = DeltaW; d *= eps; Wp += d;
        const REAL phip = eval_phi(Wp);
        phip0 = (phip - phi0) / eps;
        if (phip0 >= (REAL)0) phip0 = -(REAL)std::max((REAL)1e-16, phi0); // força descida
    }

    // se ΔW vazio, não anda
    if (DeltaW.Rows()==0 || DeltaW.Cols()==0) { NextW=Wn; this->LoadSolution(backup); return (REAL)0.0; }

    REAL alpha = (REAL)1.0;
    REAL best_phi = phi0, best_alpha = (REAL)0.0;
    TPZFMatrix<STATE> bestW = Wn;

    for (int k = 0; k < niter && alpha >= amin; ++k) {
        TPZFMatrix<STATE> trial = Wn, d = DeltaW; d *= alpha; trial += d;

        const REAL phi_a = eval_phi(trial);
        if (phi_a < best_phi) { best_phi = phi_a; best_alpha = alpha; bestW = trial; }

        // Armijo: φ(α) ≤ φ(0) + c1 α φ'(0)
        if (phi_a <= (REAL)(phi0 + c1 * alpha * phip0)) {
            NextW = trial;
            this->LoadSolution(backup);
            return alpha;
        }

        // --- interpolação quadrática protegida ---
        // modelo: φ(α) ≈ φ0 + φ'0 α + c α^2  =>  c = (φ(α) - φ0 - φ'0 α)/α^2
        REAL denom = (phi_a - phi0 - phip0*alpha);
        REAL a_quad;
        if (denom <= (REAL)0) {
            a_quad = alpha * shrink_hi; // guarda-corpo: cai para meia-passada
        } else {
            a_quad = -(phip0) * alpha * alpha / ( (REAL)2.0 * denom );
            // protege dentro de [0.1α, 0.5α]
            a_quad = std::max(shrink_lo*alpha, std::min(shrink_hi*alpha, a_quad));
        }
        alpha = a_quad;
    }

    // fallback: aceita o melhor que reduziu φ, senão não anda
    if (best_phi < phi0) { NextW = bestW; this->LoadSolution(backup); return best_alpha; }
    NextW = Wn; this->LoadSolution(backup); return (REAL)0.0;
}

// Armijo (backtracking) line search
// tol  -> alpha mínimo permitido (ex.: 1e-6)
// niter-> máximo de backtracks
// Armijo (backtracking) – versão que COMPILA com LoadSolution(TPZFMatrix<STATE>&)
REAL TPZElastoPlasticAnalysis::ArmijoLineSearch(const TPZFMatrix<STATE>& Wn,
                                          TPZFMatrix<STATE> DeltaW,
                                          TPZFMatrix<STATE>& NextW,
                                          REAL tol, int niter)
{
    // parâmetros do Armijo
    const REAL c    = (REAL)1e-4;  // decréscimo suficiente
    const REAL rho  = (REAL)0.5;   // redução do passo
    const REAL amin = std::max<REAL>(tol, (REAL)1e-8);

    // ---- utilitário: avalia ||R|| em um W qualquer ----
    auto eval_phi = [&](const TPZFMatrix<STATE>& W)->REAL {
        TPZFMatrix<STATE> tmp = W;        // cópia NÃO-CONST para LoadSolution
        this->LoadSolution(tmp);
        this->AssembleResidual();         // monta apenas o resíduo
        return Norm(this->fRhs);          // ||R||_2
    };

    // backup do estado para restaurar no retorno
    TPZFMatrix<STATE> backup = this->fSolution;

    // φ(0) em Wn
    TPZFMatrix<STATE> W0 = Wn;            // cópia NÃO-CONST
    const REAL phi0 = eval_phi(W0);

    // passo inicial
    REAL alpha = (REAL)1.0;

    // melhor já visto (fallback)
    REAL best_phi   = phi0;
    REAL best_alpha = (REAL)0.0;
    TPZFMatrix<STATE> bestW = Wn;

    // se ΔW vazio, não anda
    if (DeltaW.Rows()==0 || DeltaW.Cols()==0) {
        NextW = Wn;
        this->LoadSolution(backup);
        return (REAL)0.0;
    }

    // backtracking
    for (int k = 0; k < niter && alpha >= amin; ++k) {
        TPZFMatrix<STATE> trial = Wn;     // cópia NÃO-CONST
        TPZFMatrix<STATE> d     = DeltaW; // por valor → seguro
        d *= alpha;
        trial += d;

        const REAL phi_a = eval_phi(trial);
        if (phi_a < best_phi) { best_phi = phi_a; best_alpha = alpha; bestW = trial; }

        // critério de Armijo: phi(a) <= (1 - c a) phi(0)
        if (phi_a <= (REAL)((1.0 - c*alpha) * phi0)) {
            NextW = trial;
            this->LoadSolution(backup);   // restaura estado interno
            return alpha;
        }
        alpha *= rho;
    }

    // fallback: aceita o melhor que reduziu ||R||, senão não anda
    if (best_phi < phi0) {
        NextW = bestW;
        this->LoadSolution(backup);
        return best_alpha;
    } else {
        NextW = Wn;
        this->LoadSolution(backup);
        return (REAL)0.0;
    }
}

/*
REAL TPZElastoPlasticAnalysis::LineSearch(const TPZFMatrix<STATE>& Wn,
                                          TPZFMatrix<STATE> DeltaW,
                                          TPZFMatrix<STATE>& NextW,
                                          REAL tol, int niter)
{
    // ------- CONFIG DE VERBOSE -------
    constexpr bool kVerbose = false; // mude para false para silenciar
    auto V = [&](auto&&... xs){ if(kVerbose){ (std::cout << ... << xs); } };

    // ---- Intervalo escalar α em [0,1] ----
    REAL A = (REAL)0.0;
    REAL B = (REAL)1.0;

    // método áureo
    constexpr REAL phi  = (REAL)0.6180339887498949;  // (sqrt(5)-1)/2
    constexpr REAL phi2 = (REAL)1.0 - phi;           // ~0.381966...

    int res_evals = 0;

    // utilitário: avalia ||R(Wn + α ΔW)|| sem “sujar” estado final
    auto EvalResidualNorm = [&](REAL alpha) -> REAL {
        TPZFMatrix<STATE> trial = Wn;
        TPZFMatrix<STATE> d     = DeltaW;  // cópia local
        d *= alpha;
        trial += d;

        // salva solução corrente para restaurar depois
        TPZFMatrix<STATE> backup = fSolution;

        // IMPORTANTE: assumimos que AssembleResidual NÃO atualiza memória
        this->LoadSolution(trial);
        this->AssembleResidual();
        REAL val = Norm(fRhs);

        // restaura solução original
        this->LoadSolution(backup);
        res_evals++;
        return val;
    };

    // (opcional) medir extremos para debug
    REAL f0 = EvalResidualNorm((REAL)0.0);
    REAL f1 = EvalResidualNorm((REAL)1.0);

    // prepara pontos internos
    REAL L = A + phi2*(B - A);
    REAL M = A + phi *(B - A);

    REAL fL = EvalResidualNorm(L);
    REAL fM = EvalResidualNorm(M);

    int  it    = 0;
    REAL width = B - A;

    if(kVerbose){
        std::cout << "\n[LineSearch] -- início --\n";
        std::cout << "  tol(alpha)=" << tol << "  niter_max=" << niter << "\n";
        std::cout << "  f(0)=" << f0 << "  f(1)=" << f1 << "\n";
        std::cout << "  A="<<A<<"  B="<<B
        << "  L="<<L<<" (fL="<<fL<<")"
        << "  M="<<M<<" (fM="<<fM<<")\n";
    }

    // loop do áureo
    while (it < niter && width > tol) {
        // decide qual metade descartar
        if (fL > fM) {
            // descarta [A, L]
            A  = L;
            L  = M;
            fL = fM;
            M  = A + phi*(B - A);
            fM = EvalResidualNorm(M);
        } else {
            // descarta [M, B]
            B  = M;
            M  = L;
            fM = fL;
            L  = A + phi2*(B - A);
            fL = EvalResidualNorm(L);
        }

        width = B - A;
        if(kVerbose){
            std::cout << "  [it " << it
            << "] A="<<A<<"  B="<<B
            << "  L="<<L<<" fL="<<fL
            << "  M="<<M<<" fM="<<fM
            << "  width="<<width << "\n";
        }
        it++;
    }

    // motivo de parada (debug)
    if(kVerbose){
        if(width <= tol) std::cout << "  -> stop: width<=tol ("<<width<<" <= "<<tol<<")\n";
        else if(it >= niter) std::cout << "  -> stop: iters max ("<<it<<")\n";
    }

    // α ótimo aproximado = ponto médio do intervalo final
    REAL ALPHA = (REAL)0.5*(A + B);

    // constrói NextW = Wn + α ΔW
    NextW = Wn;
    DeltaW *= ALPHA;
    NextW += DeltaW;

    if(kVerbose){
        std::cout << "  ALPHA=" << ALPHA
        << "  iters=" << it
        << "  residual_evals=" << res_evals
        << "\n[LineSearch] -- fim --\n";
    }

    return ALPHA;
}*/

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
