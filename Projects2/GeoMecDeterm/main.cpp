// ============================
//  SINGLE FILE - COPY & PASTE
// ============================

#include <iostream>
#include <fstream>
#include <set>
#include <vector>
#include <limits>
#include <algorithm>
#include <cmath>
#include <string>

// ---------- NeoPZ core ----------
#include "pzgmesh.h"
#include "pzcmesh.h"
#include "pzvec.h"
#include "pzmanvector.h"
#include "pzgeoel.h"

#include "pzcompel.h"
#include "pzconnect.h"
#include "pzblock.h"
#include "pzintel.h"
#include "pzgeoelside.h"

// ---------- NeoPZ geo + refs ----------
#include "pzgeotriangle.h"

#include "TPZVTKGeoMesh.h"

// ---------- NeoPZ material elastoplástico ----------
#include "TPZElasticResponse.h"
#include "TPZElastoPlasticMem.h"
#include "TPZPlasticStepPV.h"
#include "TPZYCMohrCoulombPV.h"
#include "TPZMatElastoPlastic2D.h"

// ---------- NeoPZ análise ----------

#include "pzstepsolver.h"
#include <iostream>
#include <fstream>
#include <thread>
#include <vector>
#include "TPZFileStream.h"
#include <TPZBFileStream.h>
#include "SlopeAnalysis.h"
#include <fstream>
#include <iostream>
#include <fstream>
#include <thread>
#include <vector>
#include <mutex>
#include <iostream>
#include <fstream>
#include <unistd.h>
#include <sys/wait.h>
#include <memory>
#include "SlopeAnalysis.h"
#include "pznonlinanalysis.h"
#include "TPZEigenSolver.h"
#include "TPZKrylovEigenSolver.h"
#include "TPZLapackEigenSolver.h" // ou outro solver concreto
// ------------------------------------------------------------
// typedefs de material (compatível com seu projeto)
// ------------------------------------------------------------
typedef TPZPlasticStepPV<TPZYCMohrCoulombPV, TPZElasticResponse> TPlasticMC;
typedef TPZMatElastoPlastic2D<TPlasticMC, TPZElastoPlasticMem>   plasticmat;


/**
 * @brief Cria a malha geométrica 2D triangular do problema de talude/rampa.
 *
 * Elementos de contorno recebem os IDs:
 * -1 (base), -2 (direita), -3 (topo direito), -4 (topo esquerdo), -5 (esquerda), -6 (rampa).
 * Ao final, grava um VTK ("gmeshtri.vtk") com a malha geométrica.
 *
 * @param[in] ref Número de níveis de refino geométrico uniforme aplicados após a criação (>=0).
 * @return Ponteiro para TPZGeoMesh recém-criado (propriedade do chamador).
 */
TPZGeoMesh*  TriGMesh(int ref);

/**
 * @brief Cria e configura o material elastoplástico 2D (Mohr–Coulomb + elasticidade linear).
 *
 * Define coesão e ângulo de atrito (em radianos), estado plano (plane strain por padrão),
 * fator de carga 1.0 e vetor de corpo.
 *
 * @param[in] young    Módulo de Young.
 * @param[in] poisson  Coeficiente de Poisson.
 * @param[in] coes     Coesão (mesmas unidades de tensão do material).
 * @param[in] atrito   Ângulo de atrito (radianos).
 * @param[in] bodyforce Força de corpo (tamanho 3; usa-se tipicamente a componente y).
 * @param[in] planestrain 1=plane strain (padrão), 0=plane stress (se suportado pelo material).
 * @param[in] matid    ID do material (default 1).
 * @return Ponteiro para o material criado (propriedade do chamador).
 */
plasticmat*  CreateMaterial(REAL young, REAL poisson, REAL coes, REAL atrito,
                            TPZManVector<REAL,3> bodyforce,
                            int planestrain = 1, int matid = 1);

/**
 * @brief Cria a malha computacional H¹ contínua com memória, insere o material e as CCs.
 *
 * CCs impostas por IDs: -1 (fixa u_x e u_y), -2 (fixa u_x), -5 (fixa u_x).
 * Executa AutoBuild, ajusta elementos de contorno e limpa nós desconectados.
 *
 * @param[in] gmesh  Malha geométrica.
 * @param[in] pOrder Ordem polinomial preferencial.
 * @param[in] mat    Material principal (já configurado).
 * @return Ponteiro para TPZCompMesh criada (propriedade do chamador).
 */
TPZCompMesh* CreateCMesh(TPZGeoMesh* gmesh, int pOrder, plasticmat* mat);

/**
 * @brief Inicializa/repovoa a memória constitutiva em todos os pontos de integração.
 *
 * Grava em cada TPZElastoPlasticMem as propriedades fmatprop/fmatpropinit = {coesão, atrito, atrito}.
 * Usa SetUpdateMem(true) durante a escrita e desativa ao final.
 *
 * @param[in,out] cmesh  Malha computacional.
 * @param[in] coesion    Coesão a gravar nos pontos de integração.
 * @param[in] atrito     Ângulo de atrito (radianos) a gravar nos pontos de integração.
 */
void InitializeMemory(TPZCompMesh* cmesh, REAL coesion, REAL atrito);

/**
 * @brief Calcula um indicador escalar por elemento a partir da deformação plástica equivalente.
 *
 * Para cada elemento de dimensão completa, varre os pontos de integração,
 * computa \f$J_2(\varepsilon^p)\f$, usa \f$\sqrt{J_2}\f$ como medida equivalente e
 * guarda o **máximo** por elemento.
 * Publica os valores na coluna 0 de `cmesh->ElementSolution()`.
 *
 * @param[in]  cmesh                Malha computacional.
 * @param[out] fPlasticDeformSqJ2   Vetor (size = NElements) com o indicador por elemento.
 */
void ComputeElementDeformation(TPZCompMesh* cmesh, TPZVec<REAL>& fPlasticDeformSqJ2);

/**
 * @brief Divide (h-refina) elementos cujo indicador (coluna 0 de ElementSolution) excede um limiar.
 *
 * Mantém a ordem polinomial nos filhos, faz balanceamento para evitar diferença de nível > 1,
 * ajusta elementos de contorno, limpa nós desconectados, re-inicializa blocos/solução e
 * reseta a memória dos materiais com memória.
 *
 * @param[in,out] cmesh          Malha computacional.
 * @param[in]     refineaboveval Limiar do indicador para refinar.
 * @param[out]    out_newels     Conjunto com os índices dos novos elementos criados.
 */
void DivideElementsAbove(TPZCompMesh* cmesh, REAL refineaboveval, std::set<int64_t>& out_newels);

/**
 * @brief Passo de h-refino dirigido por indicador (wrapper prático).
 *
 * Recalcula o indicador, chama DivideElementsAbove() e, se houve novos elementos,
 * chama InitializeMemory() para reconstituir o estado constitutivo.
 *
 * @param[in,out] cmesh           Malha computacional.
 * @param[in]     coes            Coesão (para InitializeMemory).
 * @param[in]     atrito          Ângulo de atrito (rad) (para InitializeMemory).
 * @param[in]     refineAboveVal  Limiar do indicador para refinar.
 * @return true se houve novos elementos; false caso contrário.
 */
bool Hrefine(TPZCompMesh* cmesh, REAL coes, REAL atrito, REAL refineAboveVal);

/**
 * @brief Retorna o deslocamento vertical (u_y) no nó geométrico de coordenadas (x,y).
 *
 * Procura o elemento contendo (x,y), identifica o nó local com essas coordenadas e
 * lê o 2º DOF do connect correspondente (assumindo [u_x, u_y] por nó).
 *
 * @param[in] cmesh Malha computacional.
 * @param[in] x     Coordenada x do nó.
 * @param[in] y     Coordenada y do nó.
 * @return Valor de u_y no nó; 0.0 se não encontrado.
 */
REAL UyAtNode(TPZCompMesh* cmesh, REAL x, REAL y);

/**
 * @brief Resolve o problema não-linear para um fator de carga (FS) e aceita a solução se convergir.
 *
 * Reinicializa a memória constitutiva e zera a solução antes de resolver.
 * Usa matriz Skyline + solver direto LDL^T e `IterativeProcess` com line search habilitado.
 *
 * @param[in,out] cmesh      Malha computacional.
 * @param[in]     coes       Coesão (para InitializeMemory).
 * @param[in]     atrito     Ângulo de atrito (rad) (para InitializeMemory).
 * @param[in]     factor     Fator de carga (FS) a aplicar.
 * @param[out]    iters_out  Número de iterações executadas.
 * @return true se convergiu e `AcceptSolution()` foi chamado; false caso contrário.
 */
bool RunAndAccept(TPZCompMesh* cmesh,
                  REAL coes, REAL atrito, REAL factor, int& iters_out);

/**
 * @brief Varre FS em rampa (fs_start:fs_step:fs_max) até falhar; refina por indicador, repete (máx. max_refines).
 *
 * Guarda o último FS convergente em cada ciclo; no último ciclo, retorna o valor de `fs_max` atingido.
 * Útil como pré-processo de refinamento antes de uma busca mais fina do FS crítico.
 *
 * @param[in,out] cmesh          Malha computacional.
 * @param[in]     coes           Coesão (para InitializeMemory).
 * @param[in]     atrito         Ângulo de atrito (rad).
 * @param[in]     refineAboveVal Limiar do indicador para refinar.
 * @param[in]     max_refines    Máximo de ciclos de h-refino.
 * @param[in]     fs_start       FS inicial da rampa.
 * @param[in]     fs_step        Incremento de FS por passo.
 * @param[in]     fs_max         Teto de FS da rampa.
 * @return Último `fs_max` alcançado no ciclo final (ver logs para o último FS convergente).
 */
REAL AutoRefine(TPZCompMesh* cmesh, REAL coes, REAL atrito, REAL refineAboveVal,
                int max_refines, REAL fs_start, REAL fs_step, REAL fs_max);

/**
 * @brief Bisseção do FS crítico com robustez (garante piso que converge e teto que falha).
 *
 * Ajusta o bracket multiplicativamente até obter lo=converge e hi=falha, então
 * itera bisseção até que a largura relativa \f$(hi-lo)/(\tfrac{hi+lo}{2})\f$ seja <= tol_fs_rel.
 *
 * @param[in,out] cmesh    Malha computacional.
 * @param[in]     coes     Coesão (para InitializeMemory dentro de RunAndAccept).
 * @param[in]     atrito   Ângulo de atrito (rad).
 * @param[in]     lo       Palpite inicial para o piso (irá convergir após ajuste).
 * @param[in]     hi       Palpite inicial para o teto (irá falhar após ajuste).
 * @param[in]     tol_fs_rel Tolerância relativa da bisseção (ex.: 1e-2 para ~1%).
 * @param[in]     max_bis  Máximo de iterações de bisseção.
 * @return Estimativa do FS* na malha corrente (retorna o piso final `lo`).
 */
REAL BisectionFS(TPZCompMesh* cmesh, REAL coes, REAL atrito,
                 REAL lo, REAL hi, REAL tol_fs_rel, int max_bis);

/**
 * @brief Driver de solução com ciclos de h-refino + bisseção do FS.
 *
 * Para cada ciclo, calcula FS* por BisectionFS, imprime variação |FS−FS_old|,
 * aplica Hrefine() e repete algumas vezes, encolhendo o bracket ao redor do FS obtido.
 * (Veja logs; quebras estão comentadas para facilitar experimentos.)
 *
 * @param[in,out] cmesh Malha computacional.
 * @param[in]     coes  Coesão (para InitializeMemory).
 * @param[in]     atrito Ângulo de atrito (rad).
 */
void Solve(TPZCompMesh* cmesh, REAL coes, REAL atrito);

// ============================================================
// main
// ============================================================
int main()
{
        // 1) Malha geométrica
        int ref = 3;
        TPZGeoMesh* gmesh = TriGMesh(ref);

        // 2) Material
        REAL young   = 20000.;
        REAL poisson = 0.49;
        REAL coes    = 10.;
        REAL atrito  = 30. * M_PI / 180.;

        TPZManVector<REAL,3> bodyforce(3,0.0);
        bodyforce[1] = -20.0;

        plasticmat* mat = CreateMaterial(young, poisson, coes, atrito, bodyforce);
        int pOrder = 2;
        TPZCompMesh* cmesh = CreateCMesh(gmesh, pOrder, mat);
        mat->SetBodyForce(bodyforce);

        auto* body = dynamic_cast<plasticmat*>(cmesh->FindMaterial(1));
        if (!body) { std::cerr << "Material id=1 não encontrado.\n"; return 1; }

        int iters_out;
        REAL factor=0.8;
        REAL x=30.;
        REAL y=40.;
        for(int iload=0;iload<15;iload++)
        {
                bool converged= RunAndAccept( cmesh, coes,  atrito,  factor,  iters_out);
                REAL sol = UyAtNode(cmesh,  x,  y);
                factor+=0.1;
                cout << "uy = "<<sol <<" factor = "<< factor << " converged ="<<converged <<endl;
        }

        // Solve(cmesh, coes, atrito);
        //
        // {
        //         std::ofstream vtk1("gmeshtri_refined_preGI.vtk");
        //         TPZVTKGeoMesh::PrintGMeshVTK(cmesh->Reference(), vtk1, true);
        //         std::cout << "[VTK] gmeshtri_refined_preGI.vtk escrito.\n";
        // }


        return 0;
}

// ============================================================
// Implementações
// ============================================================
void Solve(TPZCompMesh* cmesh,REAL coes,REAL atrito)
{
        REAL FSOLD,FS;
        int neq=cmesh->NEquations();
        int neqold;
        cout << "NUMBER OF EQUATIONS  = " << neq << endl;
        for ( int iref=1; iref<=5; iref++ ) {

                cout << "# of equations  = " <<neq << " fabs(FS-FSOLD)  "  << fabs ( FS-FSOLD )  << endl;
                FSOLD=FS;
                REAL lo=0.5;
                REAL tol_fs_rel=0.01;
                int max_bis = 20;
                REAL hi=10;
                FSOLD=FS;
                FS=  BisectionFS(cmesh, coes, atrito, lo,  hi, tol_fs_rel ,  max_bis);
                neqold=neq;
                Hrefine(cmesh,coes, atrito,0.01);
                if ( fabs ( FS-FSOLD ) <0.01 ) {
                        cout << " FS-FSOLD = "<< fabs ( FS-FSOLD ) <<endl;
                        //break;
                }
                neq=cmesh->NEquations();
                if(neq==neqold)
                {
                        //break;
                }
        }
}
bool RunAndAccept(TPZCompMesh* cmesh,
                  REAL coes, REAL atrito, REAL factor, int& iters_out)
{
        auto* body = dynamic_cast<plasticmat*>(cmesh->FindMaterial(1));
        body->SetLoadFactor(factor);
        InitializeMemory(cmesh, coes, atrito);
        cmesh->Solution().Zero();

        TPZElastoPlasticAnalysis anal(cmesh, std::cout);
        TPZSkylineStructMatrix<STATE> matskl(cmesh); matskl.SetNumThreads(12);
        anal.SetStructuralMatrix(matskl);
        TPZStepSolver<STATE> step; step.SetDirect(ELDLt);
        anal.SetSolver(step);

        int iters=30;
        bool ok = anal.IterativeProcess(std::cout, (REAL)1e-2, iters, true, false, iters_out);
        if (!ok) return false;
        anal.AcceptSolution();
        //cmesh->LoadSolution(anal.CumulativeSolution());
        return true;
}

bool Hrefine(TPZCompMesh* cmesh, REAL coes, REAL atrito,REAL refineAboveVal)
{
        const int nels_before = cmesh->NElements();
        TPZVec<REAL> defel;
        ComputeElementDeformation(cmesh, defel);

        std::set<int64_t> novos;
        DivideElementsAbove(cmesh, refineAboveVal, novos);

        const int nels_after = cmesh->NElements();
        std::cout << "[PreRefine] nels: " << nels_before << " -> " << nels_after<< "  (refinados: " << (int)novos.size() << ")\n";

        if (nels_after > nels_before) {
                InitializeMemory(cmesh, coes, atrito);
                return true;
        } else {
                std::cout << "[PreRefine] sem novos refinamentos; fim.\n";
                return false;
        }
}
REAL AutoRefine(TPZCompMesh* cmesh,REAL coes, REAL atrito, REAL refineAboveVal,int  max_refines,REAL fs_start, REAL fs_step, REAL fs_max)
{

        fs_step = std::fabs(fs_step);
        if (fs_step <= (REAL)0)  fs_step = (REAL)0.05;
        if (fs_start <= (REAL)0) fs_start = (REAL)0.05;
        if (fs_max <= fs_start)  fs_max = fs_start + (REAL)10*fs_step;

        for (int ciclo = 0; ciclo < max_refines; ++ciclo) {
                std::cout << "[PreRefine] ===== ciclo " << ciclo << " =====\n";

                REAL fs = fs_start;
                REAL last_ok = (REAL)0.0;
                int  last_ok_iters = 0;

                // 1) Sobe FS até não convergir
                while (fs <= fs_max + (REAL)1e-12) {
                        int iters = 0;
                        const bool ok = RunAndAccept(cmesh, coes, atrito, fs, iters);
                        std::cout << "  [Ramp] FS=" << fs << (ok ? " (ok)" : " (falhou)")<< " | iters=" << iters << "\n";
                        if (!ok) break;
                        last_ok = fs;
                        last_ok_iters = iters;
                        fs += fs_step;
                }

                if(ciclo==max_refines-1)fs_max=fs;
                if (last_ok <= (REAL)0) {
                        std::cout << "[PreRefine] nenhum passo convergiu; encerrando.\n";
                        break;
                }
                std::cout << "[PreRefine] último FS convergente: " << last_ok<< " (iters=" << last_ok_iters << ")\n";

                bool isrefined= Hrefine(cmesh,  coes,  atrito, refineAboveVal);
                if (isrefined) {
                        InitializeMemory(cmesh, coes, atrito);
                } else {
                        std::cout << "[PreRefine] sem novos refinamentos; fim.\n";
                        break;
                }

        }

        return fs_max;

}
REAL BisectionFS(TPZCompMesh* cmesh, REAL coes, REAL atrito,
                 REAL lo, REAL hi,
                 REAL tol_fs_rel, int max_bis)
{
        auto rel_gap = [](REAL a, REAL b){
                const REAL m = (REAL)0.5*(a+b);
                return (b - a) / std::max<REAL>(m, (REAL)1e-12);
        };

        if (hi < lo) std::swap(lo, hi);

        std::cout << "\n[bisect] ===== Início da bisseção de FS =====\n"<< "[bisect] alvo: tol_rel=" << tol_fs_rel<< "  max_bis=" << max_bis<< "  lo_in=" << lo << "  hi_in=" << hi << "\n";

        // --- garantir bracket: lo converge, hi falha ---
        int it = 0, tries = 0;

        // Piso (lo) deve CONVERGIR
        std::cout << "[bisect][bracket] garantindo piso (lo) que CONVERGE...\n";
        while (!RunAndAccept(cmesh, coes, atrito, lo, it) && tries < 8) {
                std::cout << "  lo=" << lo << " -> FALHA (iters=" << it<< ")  reduzindo lo para " << (REAL)0.5*lo << "\n";
                lo *= (REAL)0.5;
                tries++;
        }
        if (tries >= 8) {
                std::cout << "[bisect][bracket][WARN] não consegui piso que converge após " << tries<< " tentativas. Prosseguindo com lo=" << lo << " (best effort).\n";
        } else {
                std::cout << "  lo=" << lo << " -> OK (iters=" << it << ")\n";
        }

        // Teto (hi) deve FALHAR
        tries = 0;
        std::cout << "[bisect][bracket] garantindo teto (hi) que FALHA...\n";
        while (RunAndAccept(cmesh, coes, atrito, hi, it) && tries < 12) {
                std::cout << "  hi=" << hi << " -> OK (iters=" << it
                << ")  aumentando hi para " << (REAL)1.5*hi << "\n";
                hi *= (REAL)1.5;
                tries++;
        }
        if (tries >= 12) {
                std::cout << "[bisect][bracket][WARN] não consegui teto que falha após " << tries
                << " tentativas. Prosseguindo com hi=" << hi << " (best effort).\n";
        } else {
                std::cout << "  hi=" << hi << " -> FALHA (iters=" << it << ")\n";
        }

        std::cout << "[bisect] bracket inicial: lo=" << lo << " (OK), hi=" << hi
        << " (FAIL)  gap_rel=" << rel_gap(lo,hi) << "\n";

        // --- bisseção ---
        int k = 0;
        while (k < max_bis) {
                const REAL gap = rel_gap(lo, hi);
                if (gap <= tol_fs_rel) {
                        std::cout << "[bisect][stop] gap_rel=" << gap
                        << " <= tol_rel=" << tol_fs_rel
                        << "  it=" << k << "\n";
                        break;
                }

                const REAL mid = (REAL)0.5*(lo + hi);
                int it_mid = 0;
                const bool ok = RunAndAccept(cmesh, coes, atrito, mid, it_mid);

                std::cout << "[bisect][it " << k << "] mid=" << mid
                << "  gap_rel=" << gap
                << "  -> " << (ok ? "OK" : "FAIL")
                << " (iters=" << it_mid << ")  ";

                if (ok) {
                        lo = mid;
                        std::cout << "novo lo=" << lo << "\n";
                } else {
                        hi = mid;
                        std::cout << "novo hi=" << hi << "\n";
                }
                k++;
        }

        const REAL fs_star = lo;
        std::cout << "[bisect] ===== FIM =====  FS*≈" << fs_star
        << "  gap_rel_final=" << rel_gap(lo,hi)
        << "  it_bis=" << k << "\n";

        return fs_star; // melhor estimativa do FS* nesta malha
}
TPZGeoMesh* TriGMesh(int ref)
{
        TPZGeoMesh* gmesh = new TPZGeoMesh();
        gmesh->SetDimension(2);

        std::vector<std::vector<double>> co = {
                /*0*/{0,0},/*1*/{10,0},/*2*/{20,0},/*3*/{30,0},/*4*/{40,0},/*5*/{50,0},/*6*/{60,0},/*7*/{70,0},
                /*8*/{0,10},/*9*/{10,10},/*10*/{20,10},/*11*/{30,10},/*12*/{40,10},/*13*/{50,10},/*14*/{60,10},/*15*/{70,10},
                /*16*/{0,20},/*17*/{10,20},/*18*/{20,20},/*19*/{30,20},/*20*/{40,20},/*21*/{50,20},/*22*/{60,20},/*23*/{70,20},
                /*24*/{0,30},/*25*/{10,30},/*26*/{20,30},/*27*/{30,30},/*28*/{40,30},/*29*/{50,30},/*30*/{60,30},/*31*/{70,30},
                /*32*/{0,40},/*33*/{10,40},/*34*/{20,40},/*35*/{30,40}
        };

        std::vector<std::vector<int>> topol = {
                /*triangles*/
                {0,1,8},{1,9,8},{1,2,9},{2,10,9},{2,3,10},{3,11,10},{3,4,11},
                {4,12,11},{4,5,12},{5,13,12},{5,6,13},{6,14,13},{6,7,14},{7,15,14},
                {8,9,16},{9,17,16},{9,10,17},{10,18,17},{10,11,18},{11,19,18},{11,12,19},
                {12,20,19},{12,13,20},{13,21,20},{13,14,21},{14,22,21},{14,15,22},{15,23,22},
                {16,17,24},{17,25,24},{17,18,25},{18,26,25},{18,19,26},{19,27,26},{19,20,27},
                {20,28,27},{20,21,28},{21,29,28},{21,22,29},{22,30,29},{22,23,30},{23,31,30},
                {24,25,32},{25,33,32},{25,26,33},{26,34,33},{26,27,34},{27,35,34},{27,28,35},
                /*lines (BCs):*/
                {0,1},{1,2},{2,3},{3,4},{4,5},{5,6},{6,7},       // -1 bottom
                {7,15},{15,23},{23,31},                          // -2 right
                {31,30},{30,29},{29,28},                         // -3 top right
                {35,34},{34,33},{33,32},                         // -4 top left
                {32,24},{24,16},{16,8},{8,0},                    // -5 left
                {28,35}                                          // -6 ramp
        };

        gmesh->NodeVec().Resize(co.size());
        TPZVec<REAL> coord(2);
        for (int i = 0; i < (int)co.size(); i++) {
                coord[0] = co[i][0];
                coord[1] = co[i][1];
                gmesh->NodeVec()[i] = TPZGeoNode(i, coord, *gmesh);
        }

        TPZVec<long> topotri(3), topoline(2);
        for (int i = 0; i < (int)topol.size(); i++) {
                if (topol[i].size() == 3) {
                        topotri[0] = topol[i][0]; topotri[1] = topol[i][1]; topotri[2] = topol[i][2];
                        new TPZGeoElRefPattern<pzgeom::TPZGeoTriangle>(i, topotri, /*matid=*/1, *gmesh);
                } else {
                        topoline[0] = topol[i][0]; topoline[1] = topol[i][1];
                        REAL x0 = co[topoline[0]][0], y0 = co[topoline[0]][1];
                        REAL xf = co[topoline[1]][0], yf = co[topoline[1]][1];
                        REAL tol = 1.e-3, L = 70, h1 = 30, h2 = 10;

                        int bcid = 0;
                        if (std::fabs(y0-0) < tol && std::fabs(yf-0) < tol)           bcid = -1; // bottom
                        else if (std::fabs(x0-L)<tol && std::fabs(xf-L)<tol)          bcid = -2; // right
                        else if (std::fabs(y0-h1)<tol && std::fabs(yf-h1)<tol)        bcid = -3; // top right
                        else if (std::fabs(y0-(h1+h2))<tol && std::fabs(yf-(h1+h2))<tol) bcid = -4; // top left
                        else if (std::fabs(x0-0)<tol && std::fabs(xf-0)<tol)          bcid = -5; // left
                        else if (std::fabs(xf-x0)>tol && std::fabs(yf-y0)>tol)        bcid = -6; // ramp
                        else {
                                std::cout << "bc element not found.\n"; DebugStop();
                        }
                        new TPZGeoElRefPattern<pzgeom::TPZGeoLinear>(i, topoline, bcid, *gmesh);
                }
        }

        gmesh->BuildConnectivity();
        for (int d = 0; d < ref; d++) {
                int nel = gmesh->NElements();
                TPZManVector<TPZGeoEl*> sub;
                for (int iel = 0; iel < nel; iel++) {
                        gmesh->ElementVec()[iel]->Divide(sub);
                }
        }

        std::ofstream vtk("gmeshtri.vtk");
        TPZVTKGeoMesh::PrintGMeshVTK(gmesh, vtk, true);
        return gmesh;
}

plasticmat* CreateMaterial(REAL young, REAL poisson, REAL coes, REAL atrito,
                           TPZManVector<REAL,3> bodyforce, int planestrain, int matid)
{
        TPZElasticResponse ER; ER.SetEngineeringData(young, poisson);

        TPlasticMC mc;
        mc.fYC.SetUp(atrito, atrito, coes, ER);
        mc.fER = ER;
        mc.SetStrengthReductionFactor(1.0);

        auto* material = new plasticmat(matid, planestrain);
        material->SetPlasticityModel(mc);
        material->SetId(matid);
        material->SetWhichLoadVector(0);
        material->SetLoadFactor(1.0);
        material->SetBodyForce(bodyforce);
        return material;
}

TPZCompMesh* CreateCMesh(TPZGeoMesh* gmesh, int pOrder, plasticmat* mat)
{
        TPZCompMesh* cmesh = new TPZCompMesh(gmesh);
        cmesh->SetDefaultOrder(pOrder);
        cmesh->SetDimModel(2);
        cmesh->InsertMaterialObject(mat);

        TPZFMatrix<STATE> val1(2,2,0.0);
        TPZManVector<STATE,2> val2(2,0.0);

        int dir = 3;
        val2[0]=1; val2[1]=1; auto* bc0 = mat->CreateBC(mat, -1, dir, val1, val2);
        val2[0]=1; val2[1]=0; auto* bc1 = mat->CreateBC(mat, -2, dir, val1, val2);
        val2[0]=1; val2[1]=0; auto* bc2 = mat->CreateBC(mat, -5, dir, val1, val2);

        cmesh->InsertMaterialObject(bc0);
        cmesh->InsertMaterialObject(bc1);
        cmesh->InsertMaterialObject(bc2);

        cmesh->SetAllCreateFunctionsContinuousWithMem();
        cmesh->AutoBuild();
        cmesh->AdjustBoundaryElements();
        cmesh->CleanUpUnconnectedNodes();
        return cmesh;
}

void InitializeMemory(TPZCompMesh* cmesh, REAL coesion, REAL atrito)
{
        auto* pMatWithMem2 = dynamic_cast<TPZMatWithMem<TPZElastoPlasticMem>*>(cmesh->MaterialVec()[1]);
        if (!pMatWithMem2) { DebugStop(); }
        pMatWithMem2->SetUpdateMem(true);

        int nels = cmesh->NElements();
        for (int iel=0; iel<nels; iel++) {
                TPZCompEl* cel = cmesh->ElementVec()[iel];
                auto* intel = dynamic_cast<TPZInterpolationSpace*>(cel);
                if (!cel || !intel) continue;
                if (dynamic_cast<TPZMatWithMem<TPZElastoPlasticMem>*>(intel->Material()) != pMatWithMem2
                        || intel->Material()->Id() < 1) continue;

                const TPZIntPoints& intpoints = intel->GetIntegrationRule();
                int nint = intpoints.NPoints();
                TPZManVector<REAL,3> q(2,0.0);

                TPZMaterialDataT<REAL> data;
                intel->InitMaterialData(data);
                data.fNeedsSol = true;

                for (int ip=0; ip<nint; ip++) {
                        REAL w; intpoints.Point(ip, q, w);
                        data.intLocPtIndex = ip;
                        intel->ComputeRequiredData(data, q);

                        int idx = data.intGlobPtIndex;
                        TPZElastoPlasticMem &mem = pMatWithMem2->MemItem(idx);
                        mem.m_elastoplastic_state.fmatprop.Resize(3);
                        mem.m_elastoplastic_state.fmatpropinit.Resize(3);

                        mem.m_elastoplastic_state.fmatpropinit[0] = coesion;
                        mem.m_elastoplastic_state.fmatpropinit[1] = atrito;
                        mem.m_elastoplastic_state.fmatpropinit[2] = atrito;

                        mem.m_elastoplastic_state.fmatprop[0] = coesion;
                        mem.m_elastoplastic_state.fmatprop[1] = atrito;
                        mem.m_elastoplastic_state.fmatprop[2] = atrito;
                }
        }
        pMatWithMem2->SetUpdateMem(false);
}

REAL UyAtNode(TPZCompMesh* cmesh, REAL x, REAL y)
{
        cmesh->LoadReferences();
        auto* gmesh = cmesh->Reference();
        if (!gmesh) return 0.0;

        TPZManVector<REAL,3> X(3,0.0); X[0]=x; X[1]=y;
        TPZManVector<REAL,3> qsi(3,0.0);
        int64_t elindex = 0;

        TPZGeoEl* gel = gmesh->FindElement(X, qsi, elindex,gmesh->Dimension());
        if (!gel || !gel->Reference()) return 0.0;

        auto* cel = gel->Reference();
        auto* intel = dynamic_cast<TPZInterpolationSpace*>(cel);
        if (!intel) return 0.0;

        int nn = gel->NNodes();
        int local = -1;
        for (int i=0;i<nn;i++){
                TPZManVector<REAL,3> co(3,0.0);
                gel->NodePtr(i)->GetCoordinates(co);
                if (std::fabs(co[0]-x) < 1e-10 && std::fabs(co[1]-y) < 1e-10){ local = i; break; }
        }
        if (local < 0) return 0.0;

        int ic = intel->ConnectIndex(local);
        if (ic < 0) return 0.0;

        TPZConnect &c = cmesh->ConnectVec()[ic];
        int64_t seq = c.SequenceNumber();
        TPZBlock &block = cmesh->Block();
        int pos = block.Position(seq);
        int ndof = block.Size(seq);

        TPZFMatrix<REAL> sol = cmesh->Solution();
        if (pos+1 >= sol.Rows() || ndof < 2) return 0.0;

        return sol(pos+1,0);
}

void ComputeElementDeformation(TPZCompMesh* cmesh, TPZVec<REAL>& fPlasticDeformSqJ2)
{
    const int64_t nelem = cmesh->NElements();
    fPlasticDeformSqJ2.resize(nelem);
    fPlasticDeformSqJ2.Fill(0.0);

    // 1 coluna para armazenar o indicador por elemento
    cmesh->ElementSolution().Redim(nelem, 1);

    for (int64_t el = 0; el < nelem; ++el) {
        TPZCompEl* cel = cmesh->ElementVec()[el];
        if (!cel) continue;
        // ignore contorno/interfaces
        if (cel->Dimension() != cmesh->Dimension()) continue;

        // material com memória
        auto* matmem = dynamic_cast<TPZMatWithMem<TPZElastoPlasticMem>*>(cel->Material());
        if (!matmem) { fPlasticDeformSqJ2[el] = 0.0; continue; }

        TPZManVector<int64_t> memindices;
        cel->GetMemoryIndices(memindices);
        if (memindices.size()==0) { fPlasticDeformSqJ2[el] = 0.0; continue; }

        REAL vmax = 0.0;
        int  npts = 0;

        for (int64_t midx : memindices) {
            if (midx < 0) continue;
            // (opcional) checagem extra:
            // if (midx >= matmem->GetMemory().size()) continue;

            const auto& mem = matmem->MemItem(midx);
            const TPZTensor<REAL> epsp = mem.m_elastoplastic_state.EpsP();

            REAL J2 = epsp.J2();
            if (J2 < (REAL)0) J2 = 0;

            const REAL eqp = std::sqrt(J2); // medida equivalente (ajuste se quiser usar sqrt(2/3)*||dev||)
            vmax = std::max(vmax, eqp);
            ++npts;
        }

        // indicador escolhido: MÁXIMO nos IPs do elemento
        // (se quiser MÉDIA, acumule e divida por npts)
        fPlasticDeformSqJ2[el] = (npts > 0) ? vmax : (REAL)0.0;
    }

    // publica na coluna 0 do ElementSolution
    cmesh->SetElementSolution(0, fPlasticDeformSqJ2);
}


void DivideElementsAbove(TPZCompMesh* cmesh, REAL refineaboveval, std::set<int64_t>& out_newels)
{
        cmesh->LoadReferences();
        std::vector<int64_t> to_divide;
        const TPZFMatrix<STATE>& elsol = cmesh->ElementSolution();

        const int64_t ne0 = cmesh->NElements();
        for (int64_t el=0; el<ne0; el++) {
                TPZCompEl* cel = cmesh->ElementVec()[el];
                if (!cel) continue;
                auto* intel = dynamic_cast<TPZInterpolationSpace*>(cel);
                if (!intel) continue;
                if (!dynamic_cast<TPZMatWithMem<TPZElastoPlasticMem>*>(cel->Material())) continue;

                if (elsol.Rows() > el && elsol.Cols() > 0) {
                        if (elsol(el,0) > refineaboveval) to_divide.push_back(el);
                }
        }

        for (auto el : to_divide) {
                TPZCompEl* cel = cmesh->ElementVec()[el]; if (!cel) continue;
                auto* intel = dynamic_cast<TPZInterpolationSpace*>(cel); if (!intel) continue;
                int p = intel->GetPreferredOrder();
                TPZStack<int64_t> sub; const int64_t idx = cel->Index();
                intel->Divide(idx, sub, /*create_boundary_elements=*/0);
                for (int i=0;i<sub.size();i++){
                        out_newels.insert(sub[i]);
                        auto* sc = cmesh->ElementVec()[sub[i]];
                        if (auto* si = dynamic_cast<TPZInterpolationSpace*>(sc)) si->SetPreferredOrder(p);
                }
        }

        // balance
        bool changed = true;
        while (changed){
                changed = false;
                std::set<int64_t> need;
                const int64_t ne = cmesh->NElements();
                for (int64_t el=0; el<ne; el++){
                        TPZCompEl* cel = cmesh->ElementVec()[el]; if (!cel) continue;
                        auto* intel = dynamic_cast<TPZInterpolationSpace*>(cel); if (!intel) continue;
                        TPZGeoEl* gel = cel->Reference(); if (!gel) continue;

                        const int ns = gel->NSides();
                        for (int s=0;s<ns;s++){
                                TPZGeoElSide gs(gel,s);
                                if (gs.Dimension() != gel->Dimension()-1) continue;
                                TPZCompElSide big = gs.LowerLevelCompElementList2(/*onlyintersect=*/1);
                                if (!big) continue;
                                TPZGeoElSide gbig(big.Reference());
                                if (gbig.Element()->Dimension() != gel->Dimension()) continue;
                                if (gel->Level() - gbig.Element()->Level() > 1){
                                        need.insert(big.Element()->Index());
                                }
                        }
                }
                for (auto el : need){
                        TPZCompEl* cel = cmesh->ElementVec()[el]; if (!cel) continue;
                        auto* intel = dynamic_cast<TPZInterpolationSpace*>(cel); if (!intel) continue;
                        if (!dynamic_cast<TPZMatWithMem<TPZElastoPlasticMem>*>(cel->Material())) continue;
                        int p = intel->GetPreferredOrder();
                        TPZStack<int64_t> sub; const int64_t idx = cel->Index();
                        intel->Divide(idx, sub, /*create_boundary_elements=*/0);
                        for (int i=0;i<sub.size();i++){
                                out_newels.insert(sub[i]);
                                auto* sc = cmesh->ElementVec()[sub[i]];
                                if (auto* si = dynamic_cast<TPZInterpolationSpace*>(sc)) si->SetPreferredOrder(p);
                        }
                        changed = true;
                }
        }

        cmesh->AdjustBoundaryElements();
        cmesh->CleanUpUnconnectedNodes();
        cmesh->InitializeBlock();
        cmesh->ExpandSolution();

        // reset memória em todos os materiais com memória
        for (auto& kv : cmesh->MaterialVec()){
                TPZMaterial* m = kv.second;
                if (auto* mm = dynamic_cast<TPZMatWithMem<TPZElastoPlasticMem>*>(m)){
                        mm->ResetMemory();
                }
        }
}





