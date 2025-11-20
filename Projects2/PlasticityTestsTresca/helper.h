#include <pzgmesh.h>
#include <pzcmesh.h>
#include <pzmanvector.h>
#include <TPZBndCondT.h>
#include <TPZLinearAnalysis.h>
#include <pzskylstrmatrix.h>
#include <pzstepsolver.h>
#include <TPZVTKGeoMesh.h>
#include <pzbuildmultiphysicsmesh.h>
#include <pzinterpolationspace.h>

#include "Elasticity/TPZElasticity2D.h"
#include "DarcyFlow/TPZDarcyFlow.h"
#include "pzlog.h"
#include <log4cxx/basicconfigurator.h>
#include "Poisson/TPZMatPoisson.h"
#include "Elasticity/TPZMatPoroElastic2DMem.h"

#include "tpzgeoelrefpattern.h"
#include "pzgeoquad.h"
#include "TPZGeoLinear.h"
#include "pzgeotriangle.h"
#include "tpzquadrilateral.h"
#include "pzgnode.h"
#include "tpzarc3d.h"
#include "pzgeopoint.h"

#include "tpzcompmeshreferred.h"
#include "TPZRefPattern.h"
#include "tpzgeoelrefpattern.h"
#include "tpzgeoblend.h"

#include "pzmultiphysicselement.h"
#include "Elasticity/TPZMatElastic2DMem.h"
#include <iostream>
#include <fstream>
#include <map>
#include <pzgraphmesh.h>
#include <TPZSSpStructMatrix.h> //symmetric sparse matrix storage
#include <pzskylstrmatrix.h> //symmetric skyline matrix storage
#include "pzsfstrmatrix.h"
#include "TPZPardisoSolver.h"
#include <Projection/TPZL2Projection.h>
#include <TPZSSpStructMatrix.h>
#include "Plasticity/TPZMatPoroElastoPlastic2DMem.h"
//#include "TPZElasticResponse.h"
#include "Plasticity/TPZElastoPlasticMem.h"
#include "Plasticity/TPZPlasticStepPV.h"
#include "Plasticity/TPZPlasticStep.h"
#include "Plasticity/TPZThermoForceA.h"
#include "Plasticity/TPZPlasticStepVoigt.h"
#include "Plasticity/TPZVonMises.h"
#include "Plasticity/TPZYCMohrCoulombPV.h"
#include "Plasticity/TPZMatElastoPlastic2D.h"

#include "Plasticity/TPZYCMohrCoulombPV2.h"

#include "TPZMaterialDataT.h"
#include "Plasticity/pzelastoplasticanalysis.h"
#include "Plasticity/TPZElasticResponse.h"
#include "pzcompelwithmem.h"
#include "pzmultiphysicscompel.h"
#include "pzgeoquad.h"          // você tem GeoQuad no print
// #include "TPZGeoTriangle.h"  // se tiver triângulos
// #include "TPZGeoLinear.h"    // se tiver elementos 1D
#include "Plasticity/TPZYCVonMisesVoigt.h"
#include "Plasticity/TPZYCTrescaVoigt.h"
#include "pzfstrmatrix.h"

#include "pzgeotetrahedra.h"
#include "pzgeopyramid.h"
#include "TPZRefPatternTools.h"
#include "pzgeoelbc.h"

#include "pzgeoquad.h"
#include "TPZGeoLinear.h"
#include "pzgeotriangle.h"
#include "tpzgeoelrefpattern.h"
#include "pzgmesh.h"
#include "pzcmesh.h"
#include "TPZVTKGeoMesh.h"
#include "tpzquadraticquad.h"

typedef TPZPlasticStepVoigt<TPZYCTrescaVoigt, TPZElasticResponse> TPlasticStepVoigtTresca;

typedef TPZPlasticStepVoigt<TPZYCVonMisesVoigt, TPZElasticResponse> TPlasticStepVoigtVM;

typedef TPZPlasticStepVoigt<TPZYCMohrCoulombPV2, TPZElasticResponse> TPlasticStepVoigtMC;

typedef TPZMatElastoPlastic2D<TPlasticStepVoigtTresca,TPZElastoPlasticMem> TMatElastoPlaticVoigtTresca;

typedef TPZMatElastoPlastic2D<TPlasticStepVoigtVM,TPZElastoPlasticMem> TMatElastoPlaticVoigtVM;

typedef TPZMatElastoPlastic2D<TPlasticStepVoigtMC,TPZElastoPlasticMem> TMatElastoPlaticMC;

#include <string>

/// Escreve nós + conectividade de um TPZGeoMesh em formato texto simples
void WriteGeoMesh(TPZGeoMesh *gmesh, const std::string &filename);

/// Lê um arquivo salvo por WriteGeoMesh e retorna um novo TPZGeoMesh
TPZGeoMesh *ReadGeoMesh(const std::string &filename);


bool HPrefine(TPZCompMesh* cmesh,REAL refineAboveVal,int porder);
void ComputeElementDeformation(TPZCompMesh* cmesh, TPZVec<REAL>& fPlasticDeformSqJ2);
void DivideElementsAbove(TPZCompMesh* cmesh, REAL refineaboveval, std::set<int64_t>& out_newels);
void PRefineElementsAbove(TPZCompMesh* cmesh,REAL refineaboveval,std::set<int64_t>& out_newels,int porder);
#include "pzgeoel.h"
#include "pzmanvector.h"

#include <fstream>
#include <iostream>
#include <sstream>

using std::cout;
using std::cerr;
using std::endl;



#include "pzgeoel.h"
#include "pzmanvector.h"

#include <fstream>
#include <iostream>
#include <sstream>

using std::cout;
using std::cerr;
using std::endl;

void WriteGeoMesh(TPZGeoMesh *gmesh, const std::string &filename)
{
        if (!gmesh) {
                cerr << "WriteGeoMesh: null gmesh\n";
                return;
        }

        std::ofstream out(filename.c_str());
        if (!out) {
                cerr << "WriteGeoMesh: could not open file " << filename << endl;
                return;
        }

        const int dim = gmesh->Dimension();
        const int64_t nnodes = gmesh->NNodes();

        // separa elementos de volume (dim) e de contorno (dim-1)
        std::vector<int64_t> volEls;
        std::vector<int64_t> bndEls;

        for (int64_t el = 0; el < gmesh->NElements(); el++) {
                TPZGeoEl *gel = gmesh->ElementVec()[el];
                if (!gel) continue;
                if (gel->HasSubElement()) continue; // só folhas

                if (gel->Dimension() == dim)
                        volEls.push_back(el);
                else if (gel->Dimension() == dim-1)
                        bndEls.push_back(el);
        }

        out << "PZGEOMESH\n";
        out << "DIM "      << dim              << "\n";
        out << "NODES "    << nnodes           << "\n";
        out << "ELEMENTS " << volEls.size()    << "\n";
        out << "BOUNDARY " << bndEls.size()    << "\n";

        // ---------------- NÓS ----------------
        out << "# id  x  y  z\n";
        for (int64_t i = 0; i < nnodes; i++) {
                TPZGeoNode &node = gmesh->NodeVec()[i];
                TPZManVector<REAL,3> x(3,0.);
                node.GetCoordinates(x);
                out << i << " " << x[0] << " " << x[1] << " " << x[2] << "\n";
        }

        // ---------------- ELEMENTOS DE VOLUME ----------------
        out << "# id  matid  typeInt  nn  nodeIds...\n";
        int64_t eid = 0;
        for (auto el : volEls) {
                TPZGeoEl *gel = gmesh->ElementVec()[el];
                const int matid   = gel->MaterialId();
                const int typeInt = static_cast<int>(gel->Type());
                const int nn      = gel->NNodes();

                out << eid++ << " " << matid << " " << typeInt << " " << nn;
                for (int in = 0; in < nn; in++) {
                        out << " " << gel->NodeIndex(in); // 0-based
                }
                out << "\n";
        }

        // ---------------- ELEMENTOS DE CONTORNO ----------------
        out << "# boundary: id  matid  typeInt  nn  nodeIds...\n";
        int64_t bid = 0;
        for (auto el : bndEls) {
                TPZGeoEl *gel = gmesh->ElementVec()[el];
                const int matid   = gel->MaterialId();
                const int typeInt = static_cast<int>(gel->Type());
                const int nn      = gel->NNodes();

                out << bid++ << " " << matid << " " << typeInt << " " << nn;
                for (int in = 0; in < nn; in++) {
                        out << " " << gel->NodeIndex(in);
                }
                out << "\n";
        }

        out << "END\n";
        std::ofstream vtk ( "WriteGeoMesh.vtk" );
        TPZVTKGeoMesh::PrintGMeshVTK ( gmesh, vtk, true );
        cout << "WriteGeoMesh: wrote " << nnodes
        << " nodes, " << volEls.size() << " volume elements and "
        << bndEls.size() << " boundary elements to " << filename << endl;
}
// -------------------- READ -----------------------------------------

static void SkipComments(std::istream &in)
{
        std::streampos pos;
        std::string line;
        while (true) {
                pos = in.tellg();
                if (!std::getline(in, line)) break;
                if (line.size() == 0) continue;
                if (line[0] == '#') continue;
                // linha útil -> volta uma posição
                in.seekg(pos);
                break;
        }
}

TPZGeoMesh *ReadGeoMesh(const std::string &filename)
{
        std::ifstream in(filename.c_str());
        if (!in) {
                cerr << "ReadGeoMesh: could not open file " << filename << endl;
                return nullptr;
        }

        std::string line, tag;

        // 1) Header
        if (!std::getline(in, line)) {
                cerr << "ReadGeoMesh: empty file\n";
                return nullptr;
        }
        if (line.find("PZGEOMESH") != 0) {
                cerr << "ReadGeoMesh: invalid header in file " << filename << endl;
                return nullptr;
        }

        int dim = 0;
        int64_t nnodes = 0, nelemVol = 0, nelemBnd = 0;

        // DIM
        std::getline(in, line);
        {
                std::istringstream iss(line);
                iss >> tag >> dim;  // "DIM dim"
        }

        // NODES
        std::getline(in, line);
        {
                std::istringstream iss(line);
                iss >> tag >> nnodes;  // "NODES nnodes"
        }

        // ELEMENTS
        std::getline(in, line);
        {
                std::istringstream iss(line);
                iss >> tag >> nelemVol;  // "ELEMENTS nelemVol"
        }

        // BOUNDARY
        std::getline(in, line);
        {
                std::istringstream iss(line);
                iss >> tag >> nelemBnd;  // "BOUNDARY nelemBnd"
        }

        TPZGeoMesh *gmesh = new TPZGeoMesh;
        gmesh->SetDimension(dim);
        gmesh->NodeVec().Resize(nnodes);

        // --------- Nós ---------
        SkipComments(in);
        for (int64_t i = 0; i < nnodes; i++) {
                if (!std::getline(in, line)) {
                        cerr << "ReadGeoMesh: unexpected EOF while reading nodes\n";
                        delete gmesh;
                        return nullptr;
                }
                if (line.size() == 0 || line[0] == '#') { i--; continue; }

                std::istringstream iss(line);
                int64_t id;
                REAL x=0., y=0., z=0.;
                iss >> id >> x >> y >> z;

                TPZGeoNode node;
                node.SetNodeId(id);
                node.SetCoord(0, x);
                node.SetCoord(1, y);
                node.SetCoord(2, z);
                gmesh->NodeVec()[id] = node;
        }

        // --------- Elementos de volume ---------
        SkipComments(in);
        for (int64_t e = 0; e < nelemVol; e++) {
                if (!std::getline(in, line)) {
                        cerr << "ReadGeoMesh: unexpected EOF while reading elements\n";
                        delete gmesh;
                        return nullptr;
                }
                if (line.size() == 0 || line[0] == '#') { e--; continue; }

                std::istringstream iss(line);
                int64_t id;
                int matid, typeInt, nn;
                iss >> id >> matid >> typeInt >> nn;

                TPZManVector<int64_t,8> nodeIdx(nn);
                for (int i = 0; i < nn; i++) {
                        iss >> nodeIdx[i];
                }

                MElementType type = static_cast<MElementType>(typeInt);
                int64_t newIndex = id;
                gmesh->CreateGeoElement(type, nodeIdx, matid, newIndex);
        }

        // --------- Elementos de contorno ---------
        SkipComments(in);
        for (int64_t e = 0; e < nelemBnd; e++) {
                if (!std::getline(in, line)) {
                        cerr << "ReadGeoMesh: unexpected EOF while reading boundary\n";
                        delete gmesh;
                        return nullptr;
                }
                if (line.size() == 0 || line[0] == '#') { e--; continue; }

                std::istringstream iss(line);
                int64_t id;
                int matid, typeInt, nn;
                iss >> id >> matid >> typeInt >> nn;

                TPZManVector<int64_t,8> nodeIdx(nn);
                for (int i = 0; i < nn; i++) {
                        iss >> nodeIdx[i];
                }

                MElementType type = static_cast<MElementType>(typeInt);
                int64_t newIndex = gmesh->NElements(); // não precisa ser 'id'
                gmesh->CreateGeoElement(type, nodeIdx, matid, newIndex);
        }

        gmesh->BuildConnectivity();
        std::ofstream vtk ( "ReadGeoMesh.vtk" );
        TPZVTKGeoMesh::PrintGMeshVTK ( gmesh, vtk, true );
        cout << "ReadGeoMesh: read " << nnodes
        << " nodes, " << nelemVol << " volume elements and "
        << nelemBnd << " boundary elements from " << filename << endl;

        return gmesh;
}

void PostProcessVariables ( TPZStack<std::string>& scal, TPZStack<std::string>& vec )
{
        scal.Push ( "StrainPlasticJ2" );
        vec.Push ( "Displacement" );
        vec.Push ( "DisplacementDoF" );
        scal.Push ( "StressXX" );
        scal.Push ( "StressYY" );
        scal.Push ( "StressJ2" );
        scal.Push ( "StressZZ" );
        //scal.Push ( "StrainPlasticXX" );
        //scal.Push ( "StrainPlasticYY" );
        scal.Push ( "StrainPlasticZZ" );
        scal.Push ( "StrainTotalZZ" );
        //scal.Push ( "StrainElasticZZ" );
        scal.Push ( "DamageVariable" );
}

void CreatePostProcessingMesh ( TPZCompMesh* cmesh,TPZPostProcAnalysis* pproc,int matid )
{
        if ( pproc->ReferenceCompMesh() != cmesh ) {
                pproc->SetCompMesh ( cmesh );
                TPZStack<std::string> scal, vec, all;
                PostProcessVariables ( scal, vec );
                for ( auto i=0; i<scal.size();  ++i ) all.Push ( scal[i] );
                for ( auto i=0; i<vec.size();   ++i ) all.Push ( vec[i] );
                TPZVec<int> matids ( 1 );
                matids[0] = matid;

                pproc->SetPostProcessVariables ( matids, all );
                TPZFStructMatrix<REAL> str ( pproc->Mesh() );
                str.SetNumThreads ( 0 );
                pproc->SetStructuralMatrix ( str );
        }
        pproc->TransferSolution();
}

void PostElastoplastic ( TPZCompMesh* cmesh,const std::string& vtkfile,int matid,int step,int dim )
{
        TPZPostProcAnalysis pproc;
        CreatePostProcessingMesh ( cmesh, &pproc, matid );
        TPZStack<std::string> scal, vec;
        PostProcessVariables ( scal, vec );
        pproc.DefineGraphMesh ( /*dim=*/dim, scal, vec, vtkfile );
        pproc.SetStep ( step );
        pproc.PostProcess ( 0 );
}
REAL UxAtNode2D ( TPZCompMesh* cmesh, REAL x, REAL y,int dir )
{
        cmesh->LoadReferences();
        auto* gmesh = cmesh->Reference();
        if ( !gmesh ) return 0.0;

        TPZManVector<REAL,3> X= {x,y,0.};
        TPZManVector<REAL,3> qsi{0.,0., 0.};
        int64_t elindex = 0;

        TPZGeoEl* gel = gmesh->FindElement ( X, qsi, elindex,2 );
        if ( !gel || !gel->Reference() ) return 0.0;

        auto* cel = gel->Reference();
        auto* intel = dynamic_cast<TPZInterpolationSpace*> ( cel );
        if ( !intel ) return 0.0;

        int nn = gel->NNodes();
        int local;
        for ( int i=0; i<nn; i++ ) {
                TPZManVector<REAL,3> co ( 3,0.0 );
                gel->NodePtr ( i )->GetCoordinates ( co );
                if ( std::fabs ( co[0]-x ) < 1e-10 && std::fabs ( co[1]-y ) < 1e-10 ) {
                        local = i;
                        break;
                }
        }
        if ( local < 0 ) return 0.0;

        int ic = intel->ConnectIndex ( local );
        if ( ic < 0 ) return 0.0;

        TPZConnect &c = cmesh->ConnectVec() [ic];
        int64_t seq = c.SequenceNumber();
        TPZBlock &block = cmesh->Block();
        int pos = block.Position ( seq );
        int ndof = block.Size ( seq );

        TPZFMatrix<REAL> sol = cmesh->Solution();


        return sol ( pos+dir,0 );
}
REAL UyAtNode3D ( TPZCompMesh* cmesh, REAL x, REAL y,REAL z )
{
        cmesh->LoadReferences();
        auto* gmesh = cmesh->Reference();
        if ( !gmesh ) return 0.0;

        TPZManVector<REAL,3> X ( 3,0.0 );
        X[0]=x;
        X[1]=y, X[2]=z;
        TPZManVector<REAL,3> qsi ( 3,0.0 );
        int64_t elindex = 0;

        TPZGeoEl* gel = gmesh->FindElement ( X, qsi, elindex,gmesh->Dimension() );
        if ( !gel || !gel->Reference() ) return 0.0;

        auto* cel = gel->Reference();
        auto* intel = dynamic_cast<TPZInterpolationSpace*> ( cel );
        if ( !intel ) return 0.0;

        int nn = gel->NNodes();
        int local = -1;
        for ( int i=0; i<nn; i++ ) {
                TPZManVector<REAL,3> co ( 3,0.0 );
                gel->NodePtr ( i )->GetCoordinates ( co );
                if ( std::fabs ( co[0]-x ) < 1e-10 && std::fabs ( co[1]-y ) < 1e-10&& std::fabs ( co[2]-z ) < 1e-10 ) {
                        local = i;
                        break;
                }
        }
        if ( local < 0 ) return 0.0;

        int ic = intel->ConnectIndex ( local );
        if ( ic < 0 ) return 0.0;

        TPZConnect &c = cmesh->ConnectVec() [ic];
        int64_t seq = c.SequenceNumber();
        TPZBlock &block = cmesh->Block();
        int pos = block.Position ( seq );
        int ndof = block.Size ( seq );

        TPZFMatrix<REAL> sol = cmesh->Solution();
        if ( pos+2 >= sol.Rows() || ndof < 3 ) return 0.0;

        return sol ( pos+2,0 );
}

// ===== utilidades simples =====
static inline STATE Dot ( const TPZFMatrix<STATE>& a, const TPZFMatrix<STATE>& b )
{
        TPZFMatrix<STATE> at,temp;
        a.Transpose ( &at );
        at.Multiply ( b,temp );
        if ( temp.Rows() >1 ) DebugStop();
        STATE s=temp ( 0,0 );
        return s;
}

// Eq. (4.123) – passo preditor (k = 1)
static STATE compute_dlambda0_riks ( const TPZFMatrix<STATE>& dwb,
                                     const TPZFMatrix<STATE>& dw,
                                     STATE L )
{
        const STATE s    = Dot ( dw, dwb );              // Δu^T * dū
        const STATE ndwb = Norm ( dwb );
        const STATE signum = ( s > 0.0 ? -1.0 : 1.0 );   // Souza Neto 4.123
        return signum * L / ndwb;
}

// Eqs. (4.116) + (4.118) – escolha da raiz para k > 1
static STATE compute_dlambda_riks ( const TPZFMatrix<STATE>& dwb,
                                    const TPZFMatrix<STATE>& dws,
                                    const TPZFMatrix<STATE>& dw,
                                    STATE L, int& rootIdx )
{
        const STATE aa = Dot ( dwb, dwb );
        TPZFMatrix<STATE> t = dw;
        t += dws;               // t = dw + dws
        const STATE bb = 2.0 * Dot ( dwb, t );
        const STATE cc = Dot ( t,t ) - L*L;

        const STATE eps = 1e-14;
        if ( aa < eps ) {                                 // cai para linear
                rootIdx = 1;
                return ( std::fabs ( bb ) > eps ) ? ( -cc/bb ) : 0.0;
        }

        STATE disc = bb*bb - 4.0*aa*cc;
        if ( disc < 0.0 ) disc = 0.0;                     // clamp numérico
        const STATE sq = std::sqrt ( disc );

        const STATE dl1 = ( -bb - sq ) / ( 2.0*aa );      // “menor”
        const STATE dl2 = ( -bb + sq ) / ( 2.0*aa );      // “maior”

        auto score = [&] ( STATE dl )->STATE {
                TPZFMatrix<STATE> x = dw;                     // Δu^(k-1)
                TPZFMatrix<STATE> tmp = dwb;
                tmp *= dl;
                x += dws;
                x += tmp;                           // Δu^(k-1)+δu*+δλ dū
                return Dot ( x, dw );                         // maximiza (4.118)
        };
        const STATE s1 = score ( dl1 ), s2 = score ( dl2 );
        if ( s1 > s2 ) {
                rootIdx = 1;
                return dl1;
        } else         {
                rootIdx = 2;
                return dl2;
        }
}


// Correções e documentação do método de Comprimento de Arco (Riks)
// Autor da correção: copilot (@copilot)
// Observação: este ficheiro contém apenas a função corrigida e documentada.
// Tip: adapte includes e nomes de tipos conforme o seu projeto (TPZElastoPlasticAnalysis, TPZFMatrix, TPZManVector, etc).

REAL IterativeProcessArcLength ( TPZElastoPlasticAnalysis &an,
                                 int loaddir,
                                 int indexbc,
                                 int nsteps,
                                 STATE bccondval,
                                 STATE lambda0,
                                 STATE L0,
                                 std::string vtkfile )
{
        // Recupera a malha/condição de contorno
        auto cmesh = an.Mesh();
        auto* bcmat = dynamic_cast<TPZBndCondT<STATE>*> ( cmesh->FindMaterial ( indexbc ) );
        if ( !bcmat ) {
                std::cout << "[ArcLength] BndCond material não encontrado (indexbc="
                << indexbc << ")\n";
                DebugStop();
        }

        auto* matmem = dynamic_cast<TPZMatWithMem<TPZElastoPlasticMem>*>(cmesh->FindMaterial(1));
        matmem->ResetMemory();
        int dim = cmesh->Dimension();
        REAL lambda = lambda0;   // lambda corrente (pode ser alterado durante os passos)
        REAL lambdan = 10000;   // lambda aceito (para rollback seguro)
        REAL L = L0;             // comprimento de arco atual

        TPZFMatrix<STATE> u_acc = an.Solution();
        u_acc.Zero();

        const STATE load0 = bccondval;
        bcmat->Val2() [loaddir] = load0;
        an.Assemble();
        TPZFMatrix<STATE> FEXT = an.Rhs();
        bcmat->Val2() [loaddir] = 0.0;
        int step = 0;


        const int max_cuts = 12;   // limite para reduzir L no mesmo step
        int cut_count = 0;

        STATE tollambda=1.e-3;
        STATE diff=100.;
        while ( step < nsteps && diff>tollambda) {

                // Estado de trabalho parte do aceito
                TPZFMatrix<STATE> u  = u_acc;
                // NOTA: NÃO redeclarar dw aqui (shadowing). Usamos dw como incremento do passo.
                // bool de convergência e norma do resíduo
                bool okconv = false;
                REAL normR  = 1e30;
                TPZFMatrix<STATE> dw = an.Solution();
                int countwhile=1;
                dw.Zero();

                // ---------- NEWTON DE COMPRIMENTO DE ARCO ----------
                {
                        const int  maxit_inner = 20;
                        const REAL etol_inner  = 1e-4;

                        TPZFMatrix<STATE> rhs, R, dws, dwb;

                        int it = 0;
                        // Mantemos um índice de "root" para funções auxiliares de escolha de raiz,
                        // se a implementação de compute_dlambda_riks exigir.
                        int rootIdx = 0;
                        //diff=sqrt((lambda-lambdan)*(lambda-lambdan));
                        std::cout << " \n step = "<< step  << " lambda-lambdan = "<< diff << std::endl;
                        while ( it < maxit_inner && normR > etol_inner ) {

                                // monta -FINT (observação: usamos bcmat->Val2()[loaddir]=0 para evitar
                                // re-aplicar o BC na montagem do interno; o externo será lambda*FEXT)
                                bcmat->Val2() [loaddir] = 0.0;
                                an.LoadSolution ( u );
                                an.Assemble();
                                rhs = an.Rhs(); // -FINT

                                // R = lambda * FEXT + (-FINT)
                                R = FEXT * lambda;
                                R += rhs;

                                // resolve para dws (incremento devido ao passo de Newton)
                                an.Rhs() = R;
                                an.Solve();
                                dws = an.Solution();

                                // resolve para dwb (direção básica / carregamento)
                                an.Rhs() = FEXT;
                                an.Solve();
                                dwb = an.Solution();

                                // atualiza dl pelo Riks
                                REAL dl = ( it == 0 ) ? compute_dlambda0_riks ( dwb, dw, L )
                                : compute_dlambda_riks ( dwb, dws, dw, L, /*rootIdx*/ * ( int[] ) {
                                        0
                                } );

                                const REAL dl_max = 1;
                                if ( dl >  dl_max ) dl =  dl_max;
                                if ( dl < -dl_max ) dl = -dl_max;

                                TPZFMatrix<STATE> dwtot = dwb*dl + dws;

                                u      += dwtot;
                                dw     += dwtot;
                                lambda += dl;

                                // remonta resíduo para checagem
                                an.LoadSolution ( u );
                                bcmat->Val2() [loaddir] = 0.0;
                                an.Assemble();
                                rhs   = an.Rhs();          // -FINT
                                R     = FEXT * lambda;     // lambda*FEXT
                                R    += rhs;               // + (-FINT)
                                normR = Norm ( R );

                                cout << " iter = " << it
                                << "  lambda = " << lambda
                                << " dl = "     << dl
                                << "  L = "     << L
                                 << " lambda-lambdan = "<< diff
                                //<< " Norm(dw) = "  << Norm(dwtot)
                                //<< " normRFEXT = "  << normRFEXT
                                << " normR = "  << normR <<  std::endl;
                                ++it;
                                //guarda de divergência: se a norma explode, aborta o inner
                                if ( normR > 2000 ) {
                                        break;
                                }
                                countwhile++;
                        } // fim do while (Newton)

                        okconv = ( normR <= etol_inner );
                }
                // ---------- FIM NEWTON ----------

                if ( okconv ) {
                        // aceita o passo
                        an.AcceptSolution(0);
                        u_acc   = an.Solution();
                        cmesh->LoadSolution(an.CumulativeSolution());
                        // salva saída VTK (função externa no seu projeto)
                        PostElastoplastic ( cmesh, vtkfile, /*matid*/1, /*out_step*/ step, dim );
                        cmesh->LoadSolution(u);
                        // pronto para o próximo step
                        ++step;
                        STATE fac = 8.  / countwhile;
                        cut_count = 0;

                        L*=fac;


                        // atualiza “lambda aceito” para rollback seguro
                        diff=sqrt((lambda-lambdan)*(lambda-lambdan));
                        std::cout << "fac = " << fac  << " lambda-lambdan = "<< diff <<  std::endl;
                        lambdan = lambda;
                } else {
                        // rollback e corta L
                        an.LoadSolution ( u_acc );
                        lambda = lambdan;   // volta lambda aceito
                        L *= 0.5;
                        ++cut_count;

                        if ( cut_count > max_cuts || L < 1e-14 ) {
                                std::cout << "[ArcLength] Falha em convergir no step " << step
                                << " após " << cut_count << " cortes de L. Abortando ciclo. "<< std::endl;;
                                break; // evita loop infinito do ciclo
                        }

                        // tenta novamente o MESMO step com L menor
                        continue;
                }

        } // fim while step < nsteps

        an.AcceptSolution();
        return 0.0;
}
bool RunAndAccept(TPZCompMesh* cmesh,REAL factor,int matid)
{
        auto* bodymat = dynamic_cast<TMatElastoPlaticVoigtTresca*>(cmesh->FindMaterial(matid));
        auto* bcmat = dynamic_cast<TPZBndCondT<STATE>*>(cmesh->FindMaterial(matid));


        TPZManVector<REAL,3> f0;
        if(bodymat)
        {
                f0=bodymat->GetBodyForce();
                TPZManVector<REAL,3> fb=f0;
                fb[1]*=factor;
                bodymat->SetBodyForce(fb);
        }else
        {
                if(!bcmat)
                {
                        std::cout << "material de contorno nao encontrado \n";
                        DebugStop();
                }
                f0=bcmat->Val2();
                bcmat->Val2()[0] *= factor;
                bcmat->Val2()[1] *= factor;


        }

        cmesh->Solution().Zero();

        TPZElastoPlasticAnalysis anal(cmesh, std::cout,TPZElastoPlasticAnalysis::ELineSearch::Dicotomic);

        if(true)
        {
                TPZFStructMatrix<REAL> str(cmesh);
                anal.SetStructuralMatrix(str);
                TPZStepSolver<REAL> direct;
                direct.SetDirect(ELU);
                anal.SetSolver(direct);
        }else{
                TPZSkylineStructMatrix<STATE> matskl(cmesh);
                matskl.SetNumThreads(16);
                anal.SetStructuralMatrix(matskl);
                TPZStepSolver<STATE> step; step.SetDirect(ELDLt);
                anal.SetSolver(step);
        }
        bool ok = anal.NewtonRaphson();
       // int iters_out;
        //bool ok = anal.IterativeProcess(std::cout, 1.e-6,100, true, false, iters_out);
        if(bodymat)
        {
                bodymat->SetBodyForce(f0);
        }else
        {
                bcmat->Val2()[0]=f0[0];
                bcmat->Val2()[1]=f0[1];
                bcmat->Val2()[2]=f0[2];
        }
        if (!ok) return false;

        anal.AcceptSolution();
        return true;
}
// Pequena utilidade para medir quão "fechado" está o bracket
static inline REAL RelGap(REAL a, REAL b) {
        const REAL m = (REAL)0.5*(a+b);
        return (b - a) / std::max<REAL>(m, (REAL)1e-12);
}
REAL FindFS_Bisection(TPZCompMesh* cmesh,
                      REAL lo, REAL hi,
                      REAL tol_fs_rel, int max_it,
                      int verbose,int loadmatid)
{
        if (hi < lo) std::swap(lo, hi);

        if (verbose) {
                std::cout << "\n[FS-Bisection] start"
                << "  lo=" << lo << "  hi=" << hi
                << "  tol_rel=" << tol_fs_rel
                << "  max_it=" << max_it << "\n";
        }

        int k = 0;
        while (k < max_it) {
                const REAL gap = RelGap(lo, hi);
                if (gap <= tol_fs_rel) {
                        if (verbose) {
                                std::cout << "[FS-Bisection] stop: gap=" << gap
                                << " <= tol=" << tol_fs_rel
                                << "  it=" << k << "\n";
                        }
                        break;
                }

                const REAL mid = (REAL)0.5*(lo + hi);
                int it_mid = 0;
                REAL resu,resf;
                const bool ok = RunAndAccept(cmesh,  mid,loadmatid);

                if (verbose) {
                        std::cout << "[FS-Bisection][it " << k << "] "
                        << "mid=" << mid
                        << " gap=" << gap
                        << " -> " << (ok ? "OK" : "FAIL")
                        << " (iters=" << it_mid << "resu =  "<<resu << "resf =  "<<resf <<  ")\n";
                }

                if (ok) lo = mid; else hi = mid;
                ++k;
        }

        if (verbose) {
                std::cout << "[FS-Bisection] end  FS≈" << lo
                << "  gap_final=" << RelGap(lo,hi)
                << "  it=" << k << "\n";
        }
        return lo; // melhor piso convergente
}

REAL Solve(TPZCompMesh* cmesh,int loadmatid,std::string vtkfile,int ref,STATE tol_fs_rel)
{
        REAL lo=0.5;
        REAL hi=30.;
        int max_bis = 20;
        int verbose=1;
        REAL FS=1000.;
        REAL FSOLD=0.;

        int porder=cmesh->GetDefaultOrder();
        std::cout << "\n[solve] ===== porder # "<< porder <<" ===== \n";
        //porder=1;
        int iters_out;
        int maxref=ref;
        TPZStack<STATE> fsstack;
        for ( int iref=1; iref<maxref; iref++ ) {

                int neq=cmesh->NEquations();
                std::cout << "\n[solve] ===== Refinamento # "<< iref <<" ====="<<" neq = " <<neq << "\n";
                FSOLD=FS;
                FS=  FindFS_Bisection(cmesh, lo,  hi, tol_fs_rel ,  max_bis,verbose,loadmatid);

                fsstack.Push(FS);
                if(FSOLD<FS)
                {
                        cout << "FSOLD<FS  "<< "FSOLD = " << FSOLD << " FS = "<< FS<<endl;
                        cout  << "FS final = " << FSOLD<<endl;
                        REAL resu,resf;
                        RunAndAccept( cmesh,  FSOLD,loadmatid);
                        return FSOLD;

                }
                //if(iref==maxref-1)break;
                //RunAndAccept( cmesh,  FS,loadmatid,false);
                HPrefine(cmesh,tol_fs_rel,porder);
                tol_fs_rel*=2;

                //porder+=1;
        }
        std::ofstream vtk("gmeshtrirefined.vtk");
        TPZVTKGeoMesh::PrintGMeshVTK(cmesh->Reference(), vtk, true);
        cout  << "FS final = " << FS<<endl;
        REAL resu,resf;
        //RunAndAccept( cmesh,  FS,loadmatid);
        return FS;

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
                        //STATE hardening= mem.m_elastoplastic_state.m_hardening;
                        //vmax = std::max(vmax, hardening);
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

void PRefineElementsAbove(TPZCompMesh* cmesh,
                          REAL refineaboveval,
                          std::set<int64_t>& out_newels,
                          int porder)
{
        if (!cmesh) return;

        // 1) Garantir referências geométricas
        cmesh->LoadReferences();

        const TPZFMatrix<STATE>& elsol = cmesh->ElementSolution();

        // 2) Sanidade: checar se ElementSolution existe
        if (elsol.Rows() == 0 || elsol.Cols() == 0) {
                std::cerr << "[PRefineElementsAbove] ElementSolution vazio: "
                << "chame quem preenche (ex.: Analysis::ComputeElementSolution) antes.\n";
                return;
        }

        const int64_t nelem = cmesh->NElements();

        for (int64_t pos = 0; pos < nelem; ++pos) {
                TPZCompEl* cel = cmesh->ElementVec()[pos];
                if (!cel) continue;

                // só elementos com espaço de interpolação (descarta especiais / multiphysics sem interp)
                auto* intel = dynamic_cast<TPZInterpolationSpace*>(cel);
                if (!intel) continue;

                // filtra materiais: só continua se for um mat com memória elastoplástica
                auto* pMatWithMem =
                dynamic_cast<TPZMatWithMem<TPZElastoPlasticMem>*>(cel->Material());
                if (!pMatWithMem) continue;

                // 3) usar o índice REAL do elemento (pode não coincidir com 'pos' se há buracos)
                const int64_t idx = cel->Index();
                if (idx < 0 || idx >= elsol.Rows()) continue;

                // 4) pega o escalar da 1ª coluna (ajuste se sua métrica estiver em outra coluna)
                const STATE val = elsol(idx, 0);
                if (val < refineaboveval) continue;

                //cout << "porder"<<porder <<endl;
                // 5) aplica p-refine: define a ordem preferida para o elemento
                intel->SetPreferredOrder(porder);

                out_newels.insert(idx);
        }

        // 6) Rebuild de conects/estrutura e vetor solução
        cmesh->AdjustBoundaryElements();
        cmesh->CleanUpUnconnectedNodes(); // remove connects órfãos
        cmesh->ExpandSolution();          // expande Solution conforme novas ordens
        cmesh->ComputeNodElCon();         // (opcional) recomputa incidências
        cmesh->InitializeBlock();         // re-inicializa blocagem conforme connects
}
bool HPrefine(TPZCompMesh* cmesh,REAL refineAboveVal,int porder)
{
        const int nels_before = cmesh->NElements();
        TPZVec<REAL> defel;
        ComputeElementDeformation(cmesh, defel);

        std::set<int64_t> novosp;
        PRefineElementsAbove (cmesh,refineAboveVal, novosp, porder );

        std::set<int64_t> novos;
        DivideElementsAbove(cmesh, refineAboveVal, novos);

        const int nels_after = cmesh->NElements();
        std::cout << "[HRefine] nels: " << nels_before << " -> " << nels_after
        << "  (refinados: " << (int)novos.size() << ")\n";

        std::cout << "[PRefine] nels: " << (int)novosp.size() << "\n";

        if (nels_after > nels_before) {
                return true;
        } else {
                std::cout << "[PreRefine] sem novos refinamentos; fim.\n";
                return false;
        }
}


void ApplyLoad ( TPZCompMesh* cmesh,TPZManVector<REAL> factors,int loaddir,int indexbc,std::string vtkfile )
{
        int dim=cmesh->Dimension();

        cout << "dimensao da malha computacional = " << dim << endl;

        // parâmetros de controle
        int nloads = factors.size();

        TPZElastoPlasticAnalysis anal ( cmesh, std::cout,TPZElastoPlasticAnalysis::ELineSearch::Armijo );

        auto* bcmat = dynamic_cast<TPZBndCondT<STATE>*> ( cmesh->FindMaterial ( indexbc ) );

        if ( !bcmat ) {
                std::cout << "material do contorno nao encontrado"<<std::endl;
                DebugStop();

        }
        const REAL load0=bcmat->Val2() [loaddir];
        cout << " load0 = "<< load0 <<endl;
        if ( true ) {
                TPZFStructMatrix<REAL> str ( cmesh );
                anal.SetStructuralMatrix ( str );
                TPZStepSolver<REAL> direct;
                direct.SetDirect ( ELU );
                anal.SetSolver ( direct );
        } else {
                TPZSkylineStructMatrix<STATE> matskl ( cmesh );
                matskl.SetNumThreads ( 16 );
                anal.SetStructuralMatrix ( matskl );
                TPZStepSolver<STATE> step;
                step.SetDirect ( ELDLt );
                anal.SetSolver ( step );
        }

        const std::string csv_path = "loadsweep.csv";
        std::ofstream csv ( csv_path );
        csv << "step,factor,uy,iters,ok\n";
        csv << std::setprecision ( 15 ) << std::scientific;

        REAL ux=0.;
        int matid=1;
        for ( int i =0; i< nloads; i++ ) {

                //bcmat->Val2()[loaddir]=load0*factors[i];
                bcmat->Val2() [loaddir]=factors[i];
                cout << "Load step =" << i<<" factor =  "<<factors[i] <<endl;
                int iters_out;

                REAL resf,resuu;
                bool ok = anal.NewtonRaphson();
                //bool ok = anal.IterativeProcess(std::cout, 1.e-6,100, true, false, iters_out);

                ux+= UxAtNode2D ( cmesh, 100.,0,0 );
                std::cout <<"ux = "<< ux <<"\n";
                TPZFMatrix<REAL> tempsol=anal.Solution();
                bcmat->Val2() [loaddir]=0;
                anal.AcceptSolution ( );

                // cmesh->LoadSolution(anal.CumulativeSolution());
                PostElastoplastic ( cmesh,vtkfile,matid,i,dim );

                //tempsol.Zero();
                //cmesh->LoadSolution(tempsol);
                //anal.LoadSolution();
        }

}
#include <fstream>
#include <sstream>
#include <unordered_map>
#include <vector>
#include <string>
#include <algorithm>

// --- structs simples para armazenar o que lemos ---
struct NodeRec {
        double x=0,y=0,z=0;
};
struct Line2  {
        long n1, n2;
};
struct Quad4  {
        long n1, n2, n3, n4;
};
struct Tri3 {
        long n1,n2,n3;
};
static inline std::string trim ( const std::string& s )
{
        size_t a = s.find_first_not_of ( " \t\r\n" );
        size_t b = s.find_last_not_of ( " \t\r\n" );
        return ( a==std::string::npos ) ?std::string() :s.substr ( a,b-a+1 );
}
#include <unordered_map>
#include <unordered_set>
#include <cstdint>
TPZGeoMesh* ReadGiDMesh ( const std::string& filename,
                          int mat2D = 1,   // matid para elementos 2D
                          int mat1D = -1 ) // matid para elementos 1D (bordas)
{
        std::ifstream in ( filename );
        if ( !in ) {
                std::cerr << "Nao consegui abrir " << filename << "\n";
                return nullptr;
        }

        std::unordered_map<long, NodeRec> nodes;
        std::vector<Line2> lines;
        std::vector<Quad4> quads;
        std::vector<Tri3>  tris;   // <<< NOVO: triângulos

        std::string line;
        while ( std::getline ( in, line ) ) {
                line = trim ( line );
                if ( line.rfind ( "MESH", 0 ) != 0 ) continue;

                // Ex.: MESH dimension 3 ElemType Linear Nnode 2
                std::stringstream ss ( line );
                std::string tok, elemType="";
                int nnode=0;
                ss >> tok;                 // MESH
                // podemos ignorar "dimension ..."
                while ( ss >> tok ) {
                        if ( tok == "ElemType" ) {
                                ss >> elemType;
                        }
                        if ( tok == "Nnode" )    {
                                ss >> nnode;
                        }
                }

                // ---- Coordinates (pode estar vazio em blocos seguintes) ----
                while ( std::getline ( in, line ) && trim ( line ) != "Coordinates" ) { /* pula */ }
                if ( trim ( line ) == "Coordinates" ) {
                        while ( std::getline ( in, line ) ) {
                                line = trim ( line );
                                if ( line == "End Coordinates" ) break;
                                if ( line.empty() ) continue;
                                std::stringstream cs ( line );
                                long id;
                                double x,y,z;
                                if ( cs >> id >> x >> y >> z ) {
                                        // guarda ou sobrescreve (se repetiu bloco com mesmo id)
                                        nodes[id] = {x,y,z};
                                }
                        }
                }

                // ---- Elements ----
                while ( std::getline ( in, line ) && trim ( line ) != "Elements" ) { /* pula */ }
                if ( trim ( line ) != "Elements" ) break;

                while ( std::getline ( in, line ) ) {
                        line = trim ( line );
                        if ( line == "End Elements" ) break;
                        if ( line.empty() ) continue;
                        std::stringstream es ( line );
                        long eid;
                        es >> eid;

                        if ( elemType == "Linear" && nnode == 2 ) {
                                long a,b;
                                es >> a >> b;
                                lines.push_back ( {a,b} );

                        } else if ( elemType == "Quadrilateral" && nnode == 4 ) {
                                long a,b,c,d;
                                es >> a >> b >> c >> d;
                                quads.push_back ( {a,b,c,d} );

                        } else if ( elemType == "Triangle" && nnode == 3 ) { // <<< NOVO
                                long a,b,c;
                                es >> a >> b >> c;
                                tris.push_back ( {a,b,c} );

                        } else {
                                // tipos não usados neste exemplo: ignore
                        }
                }
        }

        if ( nodes.empty() ) {
                std::cerr << "Arquivo nao possui bloco Coordinates valido.\n";
                return nullptr;
        }

        // --- reindexa nós (ids do arquivo) para 0..N-1 ---
        std::vector<long> ids;
        ids.reserve ( nodes.size() );
        for ( auto &kv : nodes ) ids.push_back ( kv.first );
        std::sort ( ids.begin(), ids.end() );
        std::unordered_map<long,long> mapId2Idx;
        mapId2Idx.reserve ( ids.size() );
        for ( size_t i=0; i<ids.size(); ++i ) mapId2Idx[ids[i]] = ( long ) i;

        // --- cria TPZGeoMesh ---
        TPZGeoMesh* gmesh = new TPZGeoMesh();
        gmesh->SetDimension ( 2 );
        gmesh->NodeVec().Resize ( ids.size() );

        for ( size_t i=0; i<ids.size(); ++i ) {
                TPZVec<REAL> xc ( 3,0. );
                auto &nr = nodes[ids[i]];
                xc[0]=nr.x;
                xc[1]=nr.y;
                xc[2]=nr.z;
                gmesh->NodeVec() [i] = TPZGeoNode ( ( long ) i, xc, *gmesh );
        }

        long gelid = 0;
        TPZVec<long> topol2 ( 2 ), topol3 ( 3 ), topol4 ( 4 );

        // --- escolhe matid = 1 ou 2 de acordo com o centro do elemento ---
        auto MatIdVol = [gmesh](const TPZVec<long> &topol, int nnodes)
        {
                TPZManVector<REAL,3> X(3,0.0);
                double xc = 0.0, yc = 0.0;

                for (int i = 0; i < nnodes; ++i) {
                        gmesh->NodeVec()[topol[i]].GetCoordinates(X);
                        xc += X[0];
                        yc += X[1];
                    }
                xc /= nnodes;
                yc /= nnodes;

                const double tol = 1e-8;

                // região especial: 0 < x < 0.5  e  5 < y < 5.1
                if (xc > 0.0 + tol && xc < 0.5 - tol &&
                        yc > 5.0 + tol && yc < 5.1 - tol) {
                                return 2;                           // matid especial
                            }

                return 1;                                   // matid padrão dos volumes
        };


        // ---- 2D QUADS ----
        for (const auto &q : quads) {
                topol4[0] = mapId2Idx[q.n1];
                topol4[1] = mapId2Idx[q.n2];
                topol4[2] = mapId2Idx[q.n3];
                topol4[3] = mapId2Idx[q.n4];

                int matid = MatIdVol(topol4, 4);
                new TPZGeoElRefPattern<pzgeom::TPZGeoQuad>(gelid++, topol4, matid, *gmesh);
            }

        // ---- 2D TRIÂNGULOS (se existirem) ----
        for (const auto &t : tris) {
                topol3[0] = mapId2Idx[t.n1];
                topol3[1] = mapId2Idx[t.n2];
                topol3[2] = mapId2Idx[t.n3];

                int matid = MatIdVol(topol3, 3);
                new TPZGeoElRefPattern<pzgeom::TPZGeoTriangle>(gelid++, topol3, matid, *gmesh);
            }

        // // 2D quads
        // for ( const auto &q : quads ) {
        //         topol4[0]=mapId2Idx[q.n1];
        //         topol4[1]=mapId2Idx[q.n2];
        //         topol4[2]=mapId2Idx[q.n3];
        //         topol4[3]=mapId2Idx[q.n4];
        //         new TPZGeoElRefPattern<pzgeom::TPZGeoQuad> ( gelid++, topol4, mat2D, *gmesh );
        // }
        //
        // // 2D triângulos  <<< NOVO
        // for ( const auto &t : tris ) {
        //         topol3[0]=mapId2Idx[t.n1];
        //         topol3[1]=mapId2Idx[t.n2];
        //         topol3[2]=mapId2Idx[t.n3];
        //         new TPZGeoElRefPattern<pzgeom::TPZGeoTriangle> ( gelid++, topol3, mat2D, *gmesh );
        // }

        // chave canônica (aresta não orientada) para (a,b) com a<b
        auto edge_key = [] ( long a, long b ) -> uint64_t {
                if ( a > b ) std::swap ( a, b );
                return ( ( uint64_t ) a << 32 ) | ( uint32_t ) b;
        };

        // conta incidência de arestas em elementos 2D
        std::unordered_map<uint64_t,int> edge_count;

        // QUADs
        for ( const auto &q : quads ) {
                long a = mapId2Idx[q.n1], b = mapId2Idx[q.n2];
                long c = mapId2Idx[q.n3], d = mapId2Idx[q.n4];
                edge_count[edge_key ( a,b )]++;
                edge_count[edge_key ( b,c )]++;
                edge_count[edge_key ( c,d )]++;
                edge_count[edge_key ( d,a )]++;
        }

        // TRIs  <<< NOVO
        for ( const auto &t : tris ) {
                long a = mapId2Idx[t.n1], b = mapId2Idx[t.n2], c = mapId2Idx[t.n3];
                edge_count[edge_key ( a,b )]++;
                edge_count[edge_key ( b,c )]++;
                edge_count[edge_key ( c,a )]++;
        }

        // conjunto de arestas de fronteira = arestas que aparecem 1 vez
        std::unordered_set<uint64_t> boundary_edges;
        for ( const auto &kv : edge_count ) {
                if ( kv.second == 1 ) boundary_edges.insert ( kv.first );
        }

        // 1D lines (bordas) — criar SOMENTE se for aresta de fronteira
        for ( const auto &e : lines ) {
                long a = mapId2Idx[e.n1];
                long b = mapId2Idx[e.n2];

                // pule arestas internas (compartilhadas por 2 elementos 2D)
                if ( !boundary_edges.count ( edge_key ( a,b ) ) ) continue;

                TPZManVector<REAL,3> X0 ( 3,0. ), X1 ( 3,0. );
                gmesh->NodeVec() [a].GetCoordinates ( X0 );
                gmesh->NodeVec() [b].GetCoordinates ( X1 );

                const REAL tol = 1e-8;
                const REAL Lx  = 5.0; // x da direita
                const REAL H   = 5.0; // y do topo

                auto on = [&] ( REAL v, REAL val ) {
                        return std::abs ( v - val ) <= tol;
                };
                auto in_left_open = [&] ( REAL x ) {
                        return ( x >= 0.0 - tol ) && ( x < 0.5001 - tol );
                };

                int bcid = mat1D; // default

                if ( on ( X0[0], 0.0 ) && on ( X1[0], 0.0 ) ) {
                        bcid = -1;                                  // esquerda (x=0)
                } else if ( on ( X0[1], 0.0 ) && on ( X1[1], 0.0 ) ) {
                        bcid = -2;                                  // baixo (y=0)
                } else if ( on ( X0[0], Lx ) && on ( X1[0], Lx ) ) {
                        bcid = -3;                                  // direita (x=5)
                } else if ( on ( X0[1], H ) && on ( X1[1], H ) &&
                                in_left_open ( X0[0] ) && in_left_open ( X1[0] ) ) {
                        bcid = -4;                                  // topo com 0 ≤ x < 0.5
                }

                TPZVec<long> topol2 ( 2 );
                topol2[0]=a;
                topol2[1]=b;
                new TPZGeoElRefPattern<pzgeom::TPZGeoLinear> ( gelid++, topol2, bcid, *gmesh );
                // if(bcid==-2)
                // {
                //     new TPZGeoElRefPattern<pzgeom::TPZGeoLinear>(gelid++, topol2, 1, *gmesh);
                // }
        }

        gmesh->BuildConnectivity();
        for ( int d = 0; d < 0; d++ ) {
                int nel = gmesh->NElements();
                TPZManVector<TPZGeoEl*> sub;
                for ( int iel = 0; iel < nel; iel++ ) {
                        gmesh->ElementVec() [iel]->Divide ( sub );
                }
        }
        // opcional: exporta para conferir
        std::ofstream vtk ( "gmesh_gidD.vtk" );
        TPZVTKGeoMesh::PrintGMeshVTK ( gmesh, vtk, true );

        return gmesh;
}
void ApplyLoad2 ( TPZCompMesh* cmesh,TPZManVector<REAL> factors,int loaddir,int indexbc,std::string vtkfile )
{
        int dim=cmesh->Dimension();

        cout << "dimensao da malha computacional = " << dim << endl;

        // parâmetros de controle
        int nloads = factors.size();

        TPZElastoPlasticAnalysis anal ( cmesh, std::cout,TPZElastoPlasticAnalysis::ELineSearch::GoldenSection );

        auto* bcmat = dynamic_cast<TPZBndCondT<STATE>*> ( cmesh->FindMaterial ( indexbc ) );

        if ( !bcmat ) {
                std::cout << "material do contorno nao encontrado"<<std::endl;
                DebugStop();

        }
        const REAL load0=bcmat->Val2() [loaddir];
        cout << " load0 = "<< load0 <<endl;
        if ( true ) {
                TPZFStructMatrix<REAL> str ( cmesh );
                anal.SetStructuralMatrix ( str );
                TPZStepSolver<REAL> direct;
                direct.SetDirect ( ELU );
                anal.SetSolver ( direct );
        } else {
                TPZSkylineStructMatrix<STATE> matskl ( cmesh );
                matskl.SetNumThreads ( 16 );
                anal.SetStructuralMatrix ( matskl );
                TPZStepSolver<STATE> step;
                step.SetDirect ( ELDLt );
                anal.SetSolver ( step );
        }

        const std::string csv_path = "loadsweep.csv";
        std::ofstream csv ( csv_path );
        csv << "step,factor,uy,iters,ok\n";
        csv << std::setprecision ( 15 ) << std::scientific;

        REAL u=0.;
        int matid=1;
        REAL sumfac = 0;
        for ( int i =0; i< nloads; i++ ) {

                sumfac += factors[i];
                //bcmat->Val2()[loaddir]=load0*factors[i];
                bcmat->Val2() [loaddir]=sumfac;
                cout << "Load step =" << i<<" factor =  "<<sumfac <<endl;
                int iters_out;

                STATE tol = 1.e-3;
                TPZStack<STATE> outresF;
                TPZStack<STATE> outresU;
                // bool ok = anal.NewtonRaphson(tol, outresF, outresU);
                bool ok = anal.IterativeProcess ( std::cout, tol,100, true, false, iters_out );

                u+= UxAtNode2D(cmesh, 0.,5.,1);
                std::cout <<"u = "<< u <<std::endl;
                TPZFMatrix<REAL> tempsol=anal.Solution();
                // bcmat->Val2() [loaddir]=0;
                anal.AcceptSolution ( 0 );

                //cmesh->LoadSolution(anal.CumulativeSolution());
                // cmesh->LoadSolution(anal.CumulativeSolution());
                PostElastoplastic ( cmesh,vtkfile,matid,i,dim );
                std::ofstream vtk ( "gsimpletest.vtk" );
                TPZVTKGeoMesh::PrintGMeshVTK ( cmesh->Reference(), vtk, true );
                //tempsol.Zero();
                //cmesh->LoadSolution(tempsol);
                //anal.LoadSolution();
        }

}
struct IndexesCylinder {
        // --- parâmetros mecânicos ---
public:
        int matvolume=1;
        int bcinner=-1;
        int bcouter=-2;
        int bcbottom=-3;
        int bctop=-4;

};

TPZGeoMesh* PressurizedCylinderMesh()
{

        IndexesCylinder cylindexes;
        STATE ri=100.;//mm
        STATE re=200.;//mm
        STATE h=20.;
        STATE theta=90*M_PI/180.;
        STATE s=sin ( theta );
        STATE c=cos ( theta );
        STATE s2=sin ( theta/2. );
        STATE c2=cos ( theta/2. );

        REAL co[6][2] = {

                {ri,0},{re,0},
                {ri*c,ri*s},{re*c,re*s},
                {ri*c2,ri*s2},{re*c2,re*s2}


        };

        TPZGeoEl *elvec[1];
        TPZGeoMesh *gmesh = new TPZGeoMesh();
        gmesh->SetDimension ( 2 );
        long nnode = 6;
        long nod;
        for ( nod=0; nod<nnode; nod++ ) {
                long nodind = gmesh->NodeVec().AllocateNewElement();
                TPZVec<REAL> coord ( 2 );
                coord[0] = co[nod][0];
                coord[1] = co[nod][1];
                gmesh->NodeVec() [nodind] = TPZGeoNode ( nod,coord,*gmesh );
        }

        long id=0;
        TPZVec<long> TopolQuad ( 4 );
        TopolQuad[0] = 0;
        TopolQuad[1] = 1;
        TopolQuad[2] = 3;
        TopolQuad[3] = 2;
        new TPZGeoElRefPattern< pzgeom::TPZGeoBlend< pzgeom::TPZGeoQuad> > ( id,TopolQuad,cylindexes.matvolume,*gmesh );

        id++;
        TPZVec<long> TopolArc ( 3 );
        TopolArc[0] = 0;
        TopolArc[1] = 2;
        TopolArc[2] = 4;
        new TPZGeoElRefPattern< pzgeom::TPZArc3D > ( id,TopolArc,cylindexes.bcinner,*gmesh );

        id++;
        TopolArc[0] = 1;
        TopolArc[1] = 3;
        TopolArc[2] = 5;
        new TPZGeoElRefPattern< pzgeom::TPZArc3D > ( id,TopolArc,cylindexes.bcouter,*gmesh );

        id++;
        TPZVec<long> TopolLine ( 2 );
        TopolLine[0]=0;
        TopolLine[1]=1;
        new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > ( id,TopolLine,cylindexes.bcbottom,*gmesh );

        id++;
        TopolLine[0]=3;
        TopolLine[1]=2;
        new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > ( id,TopolLine,cylindexes.bctop,*gmesh );

        cout << "b" << endl;
        gmesh->BuildConnectivity();
        // cout << "c" << endl;
        for ( int d = 0; d<3; d++ ) {
                int nel = gmesh->NElements();
                TPZManVector<TPZGeoEl *> subels;
                for ( int iel = 0; iel<nel; iel++ ) {
                        TPZGeoEl *gel = gmesh->ElementVec() [iel];
                        gel->Divide ( subels );
                }
        }

        std::ofstream files ( "teste-blend.vtk" );
        TPZVTKGeoMesh::PrintGMeshVTK ( gmesh,files,false );
        cout << "c" << endl;
        return gmesh;
}

static TPZCompMesh* CompMeshCyl ( TPZGeoMesh* gmesh )
{

        auto *mphys = new TPZCompMesh ( gmesh );
        mphys->SetDimModel ( 2 );
        mphys->SetAllCreateFunctionsContinuousWithMem();
        mphys->SetDefaultOrder ( 2 );
        auto mat = new TMatElastoPlaticVoigtTresca ( 1 );


        TPZElasticResponse ER;
        ER.SetEngineeringData ( 210,0.3 );

        TPZYCTrescaVoigt vmyc;
        const STATE sigmaY0 = 0.240;
        const STATE Hiso    = 0.;
        vmyc.SetUp ( sigmaY0,Hiso,ER );


        TPlasticStepVoigtTresca PlasticStepVoigt;



        PlasticStepVoigt.SetPlasticCriterion ( vmyc );
        PlasticStepVoigt.SetElasticResponse ( ER );



        mat->SetPlasticityModel ( PlasticStepVoigt );
        mat->SetId ( 1 );
        mphys->InsertMaterialObject ( mat ); //0
        mat->Print ( std::cout );

        TPZFMatrix<STATE> v1 ( 3,3,0. );
        TPZManVector<STATE,3> v2 ( 3,0. );
        int dirdirichlet=3,pressure=5;

        IndexesCylinder bcindexes;
        v2[0]=0.;
        v2[1]=1.;
        mphys->InsertMaterialObject ( mat->CreateBC ( mat,bcindexes.bcbottom, dirdirichlet, v1, v2 ) );

        v2[0]=1.;
        v2[1]=0.;
        mphys->InsertMaterialObject ( mat->CreateBC ( mat,bcindexes.bctop, dirdirichlet, v1, v2 ) );

        //v2[0]=-0.19209;
        v2[0]=0;
        v2[1]=0.;
        mphys->InsertMaterialObject ( mat->CreateBC ( mat,bcindexes.bcinner, pressure, v1, v2 ) );

        mphys->AutoBuild();
        mphys->AdjustBoundaryElements();
        mphys->CleanUpUnconnectedNodes();
        //mphys->Print(cout);
        return mphys;
}


void SolveCyl()
{
        auto gmesh = PressurizedCylinderMesh();
        auto cmesh = CompMeshCyl ( gmesh );

        TPZElastoPlasticAnalysis an ( cmesh, std::cout,TPZElastoPlasticAnalysis::ELineSearch::Armijo );

        TPZFStructMatrix<STATE> str ( cmesh );
        an.SetStructuralMatrix ( str );
        TPZStepSolver<REAL> direct;
        direct.SetDirect ( ELU );
        an.SetSolver ( direct );
        int itersout;

        IndexesCylinder cylindexes;
        int loaddir=0;//direcao da pressao
        //TPZManVector<REAL,11> factors={-100.,-140.,-180.,-190.,-192.};
        TPZManVector<REAL,11> factors= {-0.1,-0.14,-0.15,-0.16,-0.165,-0.17};
        //factors*=-1;
        std::string namevtk="cylinder.vtk";
        ApplyLoad ( cmesh,factors,loaddir,cylindexes.bcinner,namevtk );
}
