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

#include "Plasticity/TPZYCVonMisesVoigt.h"

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
#include "TPZHWTools.h"
using namespace std;

typedef TPZPlasticStepVoigt<TPZYCMohrCoulombPV2, TPZElasticResponse> TPlasticStepVoigtMC;
typedef TPZMatElastoPlastic2D<TPlasticStepVoigtMC,TPZElastoPlasticMem> TMatElastoPlaticMC;

typedef TPZPlasticStepPV<TPZYCMohrCoulombPV, TPZElasticResponse> plasticmorh;
typedef TPZMatElastoPlastic2D <TPZPlasticStepPV<TPZYCMohrCoulombPV, TPZElasticResponse>, TPZElastoPlasticMem> plasticmat;
bool HPrefine(TPZCompMesh* cmesh,REAL refineAboveVal,int porder);
void PostElastoplastic(TPZCompMesh* cmesh,const std::string& vtkfile,int matid,int step,int dim);
bool RunAndAccept(TPZCompMesh* cmesh,REAL factor,int matid,bool post=false);
bool RunAndAccept(TPZCompMesh* cmesh,REAL factor,int matid,bool post)
{
    auto* bodymat = dynamic_cast<TMatElastoPlaticMC*>(cmesh->FindMaterial(matid));
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

    TPZElastoPlasticAnalysis anal(cmesh, std::cout,TPZElastoPlasticAnalysis::ELineSearch::Armijo);

    if(false)
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
    bool ok = anal.NewtonRaphson(false);
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
    if(post)
    {
        int dim=2;
        int matid=1;
        PostElastoplastic(cmesh,"post.vtk",matid,0,dim);
    }
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

void PostProcessVariables(TPZStack<std::string>& scal, TPZStack<std::string>& vec)
{
    scal.Push ( "StrainPlasticJ2" );
    vec.Push ( "Displacement" );
    scal.Push ( "EBodyForce" );
    scal.Push ( "StressXX" );
    scal.Push ( "StressYY" );
    scal.Push ( "StressZZ" );
    scal.Push ( "StrainPlasticZZ" );
    scal.Push ( "StrainTotalZZ" );
    scal.Push ( "DamageVariable" );
    scal.Push ( "EOrder" );
}



void CreatePostProcessingMesh(TPZCompMesh* cmesh,TPZPostProcAnalysis* pproc,int matid)
{
    if (pproc->ReferenceCompMesh() != cmesh) {
        pproc->SetCompMesh(cmesh);
        TPZStack<std::string> scal, vec, all;
        PostProcessVariables(scal, vec);
        for (auto i=0; i<scal.size();  ++i) all.Push(scal[i]);
        for (auto i=0; i<vec.size();   ++i) all.Push(vec[i]);
        TPZVec<int> matids(1); matids[0] = matid;
        pproc->SetPostProcessVariables(matids, all);
        TPZFStructMatrix<REAL> str(pproc->Mesh());
        str.SetNumThreads(0);
        pproc->SetStructuralMatrix(str);
    }
    pproc->TransferSolution();
}

void PostElastoplastic(TPZCompMesh* cmesh,const std::string& vtkfile,int matid,int step,int dim)
{
    TPZPostProcAnalysis pproc;
    CreatePostProcessingMesh(cmesh, &pproc, matid);
    TPZStack<std::string> scal, vec;
    PostProcessVariables(scal, vec);
    pproc.DefineGraphMesh(/*dim=*/dim, scal, vec, vtkfile);
    pproc.SetStep(step);
    pproc.PostProcess(0);
}



REAL Solve(TPZCompMesh* cmesh,int loadmatid,string vtkfile,int ref,STATE tol_fs_rel)
{
    REAL lo=0.5;
    REAL hi=30.;
    int max_bis = 20;
    int verbose=1;
    REAL FS=1000.;
    REAL FSOLD=0.;

    int porder=cmesh->GetDefaultOrder();
    porder+=1;
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
        HPrefine(cmesh,tol_fs_rel,porder+1);

        //porder+=1;
    }
    std::ofstream vtk("gmeshtrirefined.vtk");
    TPZVTKGeoMesh::PrintGMeshVTK(cmesh->Reference(), vtk, true);
    cout  << "FS final = " << FS<<endl;
    REAL resu,resf;
    //RunAndAccept( cmesh,  FS,loadmatid);
    return FS;

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
TPZCompMesh* CreateCMesh(TPZGeoMesh* gmesh, int pOrder)
{
    TPZCompMesh* cmesh = new TPZCompMesh(gmesh);
    cmesh->SetDefaultOrder(pOrder);
    cmesh->SetDimModel(2);

    STATE phi=30*M_PI/180.;
    STATE psi=phi;
    STATE c =10.;
    TPZElasticResponse ER;
    ER.SetEngineeringData(20000.,0.49);
    auto mc = TPZYCMohrCoulombPV2( phi, psi, c,ER) ;
    TPlasticStepVoigtMC PlasticStepVoigt;
    PlasticStepVoigt.SetPlasticCriterion(mc);
    PlasticStepVoigt.SetElasticResponse(ER);


    int id=1;
    int planestrain=1;
    auto* material = new TMatElastoPlaticMC(id, planestrain);
    material->SetPlasticityModel(PlasticStepVoigt);
    material->SetId(id);
    TPZManVector<REAL,3> fb ={0,-20,0};
    material->SetBodyForce0(fb);
    material->SetBodyForce(fb);

    cmesh->InsertMaterialObject(material);

    TPZFMatrix<STATE> val1(2,2,0.0);
    TPZManVector<STATE,2> val2(2,0.0);

    int dir = 3;
    val2[0]=1; val2[1]=1; auto* bc0 = material->CreateBC(material, -1, dir, val1, val2);
    val2[0]=1; val2[1]=0; auto* bc1 = material->CreateBC(material, -2, dir, val1, val2);
    val2[0]=1; val2[1]=0; auto* bc2 = material->CreateBC(material, -5, dir, val1, val2);

    cmesh->InsertMaterialObject(bc0);
    cmesh->InsertMaterialObject(bc1);
    cmesh->InsertMaterialObject(bc2);

    cmesh->SetAllCreateFunctionsContinuousWithMem();
    cmesh->AutoBuild();
    cmesh->AdjustBoundaryElements();
    cmesh->CleanUpUnconnectedNodes();
    return cmesh;
}

#include "TPZVTKGeoMesh.h"
#include "pzmanvector.h"

#include <fstream>
#include <sstream>
#include <unordered_map>
#include <vector>
#include <string>
#include <algorithm>

// --- structs simples para armazenar o que lemos ---
struct NodeRec { double x=0,y=0,z=0; };
struct Line2  { long n1, n2; };
struct Quad4  { long n1, n2, n3, n4; };
struct Tri3 {long n1,n2,n3;};
static inline std::string trim(const std::string& s){
    size_t a = s.find_first_not_of(" \t\r\n");
    size_t b = s.find_last_not_of(" \t\r\n");
    return (a==std::string::npos)?std::string():s.substr(a,b-a+1);
}
#include <unordered_map>
#include <unordered_set>
#include <cstdint>
// Lê um arquivo GiD .msh com blocos "MESH ...", "Coordinates", "Elements"
TPZGeoMesh* ReadGiDMeshSimple(const std::string& filename,
                              int mat2D = 1,   // matid para elementos 2D
                              int mat1D = -1)  // matid para elementos 1D (bordas)
{
    std::ifstream in(filename);
    if(!in){ std::cerr << "Nao consegui abrir " << filename << "\n"; return nullptr; }

    std::unordered_map<long, NodeRec> nodes;
    std::vector<Line2> lines;
    std::vector<Quad4> quads;

    std::string line;
    while (std::getline(in, line)) {
        line = trim(line);
        if (line.rfind("MESH", 0) != 0) continue;

        // Ex.: MESH dimension 3 ElemType Linear Nnode 2
        std::stringstream ss(line);
        std::string tok, elemType=""; int nnode=0;
        ss >> tok;                 // MESH
        // podemos ignorar "dimension ..."
        while (ss >> tok) {
            if (tok == "ElemType") { ss >> elemType; }
            if (tok == "Nnode")    { ss >> nnode;    }
        }

        // ---- Coordinates (pode estar vazio em blocos seguintes) ----
        while (std::getline(in, line) && trim(line) != "Coordinates") { /* pula */ }
        if (trim(line) == "Coordinates") {
            while (std::getline(in, line)) {
                line = trim(line);
                if (line == "End Coordinates") break;
                if (line.empty()) continue;
                std::stringstream cs(line);
                long id; double x,y,z;
                if (cs >> id >> x >> y >> z) {
                    // guarda ou sobrescreve (se repetiu bloco com mesmo id)
                    nodes[id] = {x,y,z};
                }
            }
        }

        // ---- Elements ----
        while (std::getline(in, line) && trim(line) != "Elements") { /* pula */ }
        if (trim(line) != "Elements") break;

        while (std::getline(in, line)) {
            line = trim(line);
            if (line == "End Elements") break;
            if (line.empty()) continue;
            std::stringstream es(line);
            long eid; es >> eid;
            if (elemType == "Linear" && nnode == 2) {
                long a,b; es >> a >> b;
                lines.push_back({a,b});
            } else if (elemType == "Quadrilateral" && nnode == 4) {
                long a,b,c,d; es >> a >> b >> c >> d;
                quads.push_back({a,b,c,d});
            } else {
                // tipos não usados neste exemplo: ignore
            }
        }
    }

    if (nodes.empty()) {
        std::cerr << "Arquivo nao possui bloco Coordinates valido.\n";
        return nullptr;
    }

    // --- reindexa nós (ids do arquivo) para 0..N-1 ---
    std::vector<long> ids; ids.reserve(nodes.size());
    for (auto &kv : nodes) ids.push_back(kv.first);
    std::sort(ids.begin(), ids.end());
    std::unordered_map<long,long> mapId2Idx; mapId2Idx.reserve(ids.size());
    for (size_t i=0;i<ids.size();++i) mapId2Idx[ids[i]] = (long)i;

    // --- cria TPZGeoMesh ---
    TPZGeoMesh* gmesh = new TPZGeoMesh();
    gmesh->SetDimension(2);
    gmesh->NodeVec().Resize(ids.size());

    for (size_t i=0;i<ids.size();++i){
        TPZVec<REAL> xc(3,0.);
        auto &nr = nodes[ids[i]];
        xc[0]=nr.x; xc[1]=nr.y; xc[2]=nr.z;
        gmesh->NodeVec()[i] = TPZGeoNode((long)i, xc, *gmesh);
    }

    long gelid = 0;
    TPZVec<long> topol2(2), topol4(4);

    // 2D quads
    for (const auto &q : quads){
        topol4[0]=mapId2Idx[q.n1];
        topol4[1]=mapId2Idx[q.n2];
        topol4[2]=mapId2Idx[q.n3];
        topol4[3]=mapId2Idx[q.n4];
        new TPZGeoElRefPattern<pzgeom::TPZGeoQuad>(gelid++, topol4, mat2D, *gmesh);
    }
    // chave canônica (aresta não orientada) para (a,b) com a<b
    auto edge_key = [](long a, long b) -> uint64_t {
        if (a > b) std::swap(a, b);
        return ( (uint64_t)a << 32 ) | (uint32_t)b;
    };

    // conta incidência de arestas em elementos 2D
    std::unordered_map<uint64_t,int> edge_count;

    // QUADs (se você tem 'quads' lidos do arquivo)
    for (const auto &q : quads) {
        long a = mapId2Idx[q.n1], b = mapId2Idx[q.n2];
        long c = mapId2Idx[q.n3], d = mapId2Idx[q.n4];
        edge_count[edge_key(a,b)]++;
        edge_count[edge_key(b,c)]++;
        edge_count[edge_key(c,d)]++;
        edge_count[edge_key(d,a)]++;
    }


    // conjunto de arestas de fronteira = arestas que aparecem 1 vez
    std::unordered_set<uint64_t> boundary_edges;
    for (const auto &kv : edge_count) {
        if (kv.second == 1) boundary_edges.insert(kv.first);
    }

    // 1D lines (bordas) — criar SOMENTE se for aresta de fronteira
    for (const auto &e : lines) {
        long a = mapId2Idx[e.n1];
        long b = mapId2Idx[e.n2];

        // pule arestas internas (compartilhadas por 2 elementos 2D)
        if (!boundary_edges.count(edge_key(a,b))) continue;

        TPZManVector<REAL,3> X0(3,0.), X1(3,0.);
        gmesh->NodeVec()[a].GetCoordinates(X0);
        gmesh->NodeVec()[b].GetCoordinates(X1);

        const REAL tol = 1e-8;
        const REAL Lx  = 5.0; // x da direita
        const REAL H   = 5.0; // y do topo

        auto on = [&](REAL v, REAL val){ return std::abs(v - val) <= tol; };
        auto in_left_open = [&](REAL x){ return (x >= 0.0 - tol) && (x < 0.5001 - tol); };

        int bcid = mat1D; // default

        if (on(X0[0], 0.0) && on(X1[0], 0.0)) {
            bcid = -1;                                  // esquerda (x=0)
        } else if (on(X0[1], 0.0) && on(X1[1], 0.0)) {
            bcid = -2;                                  // baixo (y=0)
        } else if (on(X0[0], Lx) && on(X1[0], Lx)) {
            bcid = -3;                                  // direita (x=5)
        } else if (on(X0[1], H) && on(X1[1], H) &&
            in_left_open(X0[0]) && in_left_open(X1[0])) {
            bcid = -4;                                  // topo com 0 ≤ x < 0.5
            }

            TPZVec<long> topol2(2); topol2[0]=a; topol2[1]=b;
        new TPZGeoElRefPattern<pzgeom::TPZGeoLinear>(gelid++, topol2, bcid, *gmesh);
        // if(bcid==-2)
        // {
        //     new TPZGeoElRefPattern<pzgeom::TPZGeoLinear>(gelid++, topol2, 1, *gmesh);
        // }
    }



    gmesh->BuildConnectivity();
    for (int d = 0; d < 0; d++) {
        int nel = gmesh->NElements();
        TPZManVector<TPZGeoEl*> sub;
        for (int iel = 0; iel < nel; iel++) {
            gmesh->ElementVec()[iel]->Divide(sub);
        }
    }
    // opcional: exporta para conferir
    std::ofstream vtk("gmesh_gid_simple.vtk");
    TPZVTKGeoMesh::PrintGMeshVTK(gmesh, vtk, true);

    return gmesh;
}


TPZGeoMesh* ReadGiDMesh(const std::string& filename,
                              int mat2D = 1,   // matid para elementos 2D
                              int mat1D = -1)  // matid para elementos 1D (bordas)
{
    std::ifstream in(filename);
    if(!in){
        std::cerr << "Nao consegui abrir " << filename << "\n";
        return nullptr;
    }

    std::unordered_map<long, NodeRec> nodes;
    std::vector<Line2> lines;
    std::vector<Quad4> quads;
    std::vector<Tri3>  tris;   // <<< NOVO: triângulos

    std::string line;
    while (std::getline(in, line)) {
        line = trim(line);
        if (line.rfind("MESH", 0) != 0) continue;

        // Ex.: MESH dimension 3 ElemType Linear Nnode 2
        std::stringstream ss(line);
        std::string tok, elemType=""; int nnode=0;
        ss >> tok;                 // MESH
        // podemos ignorar "dimension ..."
        while (ss >> tok) {
            if (tok == "ElemType") { ss >> elemType; }
            if (tok == "Nnode")    { ss >> nnode;    }
        }

        // ---- Coordinates (pode estar vazio em blocos seguintes) ----
        while (std::getline(in, line) && trim(line) != "Coordinates") { /* pula */ }
        if (trim(line) == "Coordinates") {
            while (std::getline(in, line)) {
                line = trim(line);
                if (line == "End Coordinates") break;
                if (line.empty()) continue;
                std::stringstream cs(line);
                long id; double x,y,z;
                if (cs >> id >> x >> y >> z) {
                    // guarda ou sobrescreve (se repetiu bloco com mesmo id)
                    nodes[id] = {x,y,z};
                }
            }
        }

        // ---- Elements ----
        while (std::getline(in, line) && trim(line) != "Elements") { /* pula */ }
        if (trim(line) != "Elements") break;

        while (std::getline(in, line)) {
            line = trim(line);
            if (line == "End Elements") break;
            if (line.empty()) continue;
            std::stringstream es(line);
            long eid; es >> eid;

            if (elemType == "Linear" && nnode == 2) {
                long a,b; es >> a >> b;
                lines.push_back({a,b});

            } else if (elemType == "Quadrilateral" && nnode == 4) {
                long a,b,c,d; es >> a >> b >> c >> d;
                quads.push_back({a,b,c,d});

            } else if (elemType == "Triangle" && nnode == 3) {   // <<< NOVO
                long a,b,c; es >> a >> b >> c;
                tris.push_back({a,b,c});

            } else {
                // tipos não usados neste exemplo: ignore
            }
        }
    }

    if (nodes.empty()) {
        std::cerr << "Arquivo nao possui bloco Coordinates valido.\n";
        return nullptr;
    }

    // --- reindexa nós (ids do arquivo) para 0..N-1 ---
    std::vector<long> ids; ids.reserve(nodes.size());
    for (auto &kv : nodes) ids.push_back(kv.first);
    std::sort(ids.begin(), ids.end());
    std::unordered_map<long,long> mapId2Idx; mapId2Idx.reserve(ids.size());
    for (size_t i=0;i<ids.size();++i) mapId2Idx[ids[i]] = (long)i;

    // --- cria TPZGeoMesh ---
    TPZGeoMesh* gmesh = new TPZGeoMesh();
    gmesh->SetDimension(2);
    gmesh->NodeVec().Resize(ids.size());

    for (size_t i=0;i<ids.size();++i){
        TPZVec<REAL> xc(3,0.);
        auto &nr = nodes[ids[i]];
        xc[0]=nr.x; xc[1]=nr.y; xc[2]=nr.z;
        gmesh->NodeVec()[i] = TPZGeoNode((long)i, xc, *gmesh);
    }

    long gelid = 0;
    TPZVec<long> topol2(2), topol3(3), topol4(4);

    // 2D quads
    for (const auto &q : quads){
        topol4[0]=mapId2Idx[q.n1];
        topol4[1]=mapId2Idx[q.n2];
        topol4[2]=mapId2Idx[q.n3];
        topol4[3]=mapId2Idx[q.n4];
        new TPZGeoElRefPattern<pzgeom::TPZGeoQuad>(gelid++, topol4, mat2D, *gmesh);
    }

    // 2D triângulos  <<< NOVO
    for (const auto &t : tris){
        topol3[0]=mapId2Idx[t.n1];
        topol3[1]=mapId2Idx[t.n2];
        topol3[2]=mapId2Idx[t.n3];
        new TPZGeoElRefPattern<pzgeom::TPZGeoTriangle>(gelid++, topol3, mat2D, *gmesh);
    }

    // chave canônica (aresta não orientada) para (a,b) com a<b
    auto edge_key = [](long a, long b) -> uint64_t {
        if (a > b) std::swap(a, b);
        return ( (uint64_t)a << 32 ) | (uint32_t)b;
    };

    // conta incidência de arestas em elementos 2D
    std::unordered_map<uint64_t,int> edge_count;

    // QUADs
    for (const auto &q : quads) {
        long a = mapId2Idx[q.n1], b = mapId2Idx[q.n2];
        long c = mapId2Idx[q.n3], d = mapId2Idx[q.n4];
        edge_count[edge_key(a,b)]++;
        edge_count[edge_key(b,c)]++;
        edge_count[edge_key(c,d)]++;
        edge_count[edge_key(d,a)]++;
    }

    // TRIs  <<< NOVO
    for (const auto &t : tris) {
        long a = mapId2Idx[t.n1], b = mapId2Idx[t.n2], c = mapId2Idx[t.n3];
        edge_count[edge_key(a,b)]++;
        edge_count[edge_key(b,c)]++;
        edge_count[edge_key(c,a)]++;
    }

    // conjunto de arestas de fronteira = arestas que aparecem 1 vez
    std::unordered_set<uint64_t> boundary_edges;
    for (const auto &kv : edge_count) {
        if (kv.second == 1) boundary_edges.insert(kv.first);
    }

    // 1D lines (bordas) — criar SOMENTE se for aresta de fronteira
    for (const auto &e : lines) {
        long a = mapId2Idx[e.n1];
        long b = mapId2Idx[e.n2];

        // pule arestas internas (compartilhadas por 2 elementos 2D)
        if (!boundary_edges.count(edge_key(a,b))) continue;

        TPZManVector<REAL,3> X0(3,0.), X1(3,0.);
        gmesh->NodeVec()[a].GetCoordinates(X0);
        gmesh->NodeVec()[b].GetCoordinates(X1);

        const REAL tol = 1e-8;
        const REAL Lx  = 5.0; // x da direita
        const REAL H   = 5.0; // y do topo

        auto on = [&](REAL v, REAL val){ return std::abs(v - val) <= tol; };
        auto in_left_open = [&](REAL x){ return (x >= 0.0 - tol) && (x < 0.5001 - tol); };

        int bcid = mat1D; // default

        if (on(X0[0], 0.0) && on(X1[0], 0.0)) {
            bcid = -1;                                  // esquerda (x=0)
        } else if (on(X0[1], 0.0) && on(X1[1], 0.0)) {
            bcid = -2;                                  // baixo (y=0)
        } else if (on(X0[0], Lx) && on(X1[0], Lx)) {
            bcid = -3;                                  // direita (x=5)
        } else if (on(X0[1], H) && on(X1[1], H) &&
                   in_left_open(X0[0]) && in_left_open(X1[0])) {
            bcid = -4;                                  // topo com 0 ≤ x < 0.5
        }

        TPZVec<long> topol2(2); topol2[0]=a; topol2[1]=b;
        new TPZGeoElRefPattern<pzgeom::TPZGeoLinear>(gelid++, topol2, bcid, *gmesh);
        // if(bcid==-2)
        // {
        //     new TPZGeoElRefPattern<pzgeom::TPZGeoLinear>(gelid++, topol2, 1, *gmesh);
        // }
    }

    gmesh->BuildConnectivity();
    for (int d = 0; d < 0; d++) {
        int nel = gmesh->NElements();
        TPZManVector<TPZGeoEl*> sub;
        for (int iel = 0; iel < nel; iel++) {
            gmesh->ElementVec()[iel]->Divide(sub);
        }
    }
    // opcional: exporta para conferir
    std::ofstream vtk("gmesh_gid_simple.vtk");
    TPZVTKGeoMesh::PrintGMeshVTK(gmesh, vtk, true);

    return gmesh;
}


TPZCompMesh* CreateCMeshFoot(TPZGeoMesh* gmesh, int pOrder)
{
    TPZCompMesh* cmesh = new TPZCompMesh(gmesh);
    cmesh->SetDefaultOrder(pOrder);
    cmesh->SetDimModel(2);

    STATE phi=20*M_PI/180.;
    STATE psi=phi;
    STATE c =490.;
    TPZElasticResponse ER;
    ER.SetEngineeringData(0.1e8,0.48);
    auto mc = TPZYCMohrCoulombPV2( phi, psi, c,ER) ;
    TPlasticStepVoigtMC PlasticStepVoigt;
    PlasticStepVoigt.SetPlasticCriterion(mc);
    PlasticStepVoigt.SetElasticResponse(ER);


    int id=1;
    int planestrain=1;
    auto* material = new TMatElastoPlaticMC(id, planestrain);
    material->SetPlasticityModel(PlasticStepVoigt);
    material->SetId(id);
    TPZManVector<REAL,3> fb ={0,0,0};
    material->SetBodyForce0(fb);
    material->SetBodyForce(fb);

    cmesh->InsertMaterialObject(material);

    TPZFMatrix<STATE> val1(2,2,0.0);
    TPZManVector<STATE,2> val2(2,0.0);

    int dir = 3;
    val2[0]=1;
    val2[1]=0;
    auto* bc0 = material->CreateBC(material, -1, dir, val1, val2);
    val2[0]=0;
    val2[1]=1;
    auto* bc1 = material->CreateBC(material, -2, dir, val1, val2);
    val2[0]=1;
    val2[1]=0;
    auto* bc2 = material->CreateBC(material, -3, dir, val1, val2);
    val2[0]=0;
    val2[1]=-490;
    auto* bc4 = material->CreateBC(material, -4, 1, val1, val2);
    // val2[0]=0;
    // val2[1]=0;
    // auto* bc4 = material->CreateBC(material, -4, 0, val1, val2);

    cmesh->InsertMaterialObject(bc0);
    cmesh->InsertMaterialObject(bc1);
    cmesh->InsertMaterialObject(bc2);
    cmesh->InsertMaterialObject(bc4);

    cmesh->SetAllCreateFunctionsContinuousWithMem();
    cmesh->AutoBuild();
    cmesh->AdjustBoundaryElements();
    cmesh->CleanUpUnconnectedNodes();
    return cmesh;
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
            STATE hardening= mem.m_elastoplastic_state.m_hardening;
            vmax = std::max(vmax, hardening);
            // REAL J2 = epsp.J2();
            // if (J2 < (REAL)0) J2 = 0;
            //
            // const REAL eqp = std::sqrt(J2); // medida equivalente (ajuste se quiser usar sqrt(2/3)*||dev||)
            // vmax = std::max(vmax, eqp);
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
REAL UyAtNode(TPZCompMesh* cmesh, REAL x, REAL y)
{
    int dim=2;
    if (cmesh->Reference()->Dimension()!=dim)
    {
        std::cout << "(cmesh->Reference()->Dimension()!=dim) = "<< (cmesh->Reference()->Dimension()!=dim) << std::endl;
        DebugStop();
    }
    cmesh->LoadReferences();
    auto* gmesh = cmesh->Reference();
    if (!gmesh) return 0.0;

    TPZManVector<REAL,3> X(3,0.0); X[0]=x; X[1]=y;
    TPZManVector<REAL,3> qsi(3,0.0);
    int64_t elindex = 0;

    TPZGeoEl* gel = gmesh->FindElement(X, qsi, elindex,dim);
    //std::cout << "elindex = "<< elindex << std::endl;
    if (!gel || !gel->Reference()) {
        std::cout << "!gel || !gel->Reference() = "<< (!gel || !gel->Reference()) << std::endl;
        DebugStop();
    }

    auto* cel = gel->Reference();
    auto* intel = dynamic_cast<TPZInterpolationSpace*>(cel);
    if (!intel)
    {
        std::cout << "!intel = "<< (!intel) << std::endl;
        DebugStop();
    }

    int nn = gel->NNodes();

    int local;
    for (int i=0;i<nn;i++){
        TPZManVector<REAL,3> co(3,0.0);
        gel->NodePtr(i)->GetCoordinates(co);
        if (std::fabs(co[0]-x) < 1e-10 && std::fabs(co[1]-y) < 1e-10){ local = i; break; }
    }
    if (local < 0)
    {
        std::cout << "local = "<<local << std::endl;
        DebugStop();
    }

    int ic = intel->ConnectIndex(local);
    if (ic < 0)
    {
        std::cout << "ic = "<<ic << std::endl;
        DebugStop();
    }

    TPZConnect &c = cmesh->ConnectVec()[ic];
    int64_t seq = c.SequenceNumber();
    TPZBlock &block = cmesh->Block();
    int pos = block.Position(seq);
    int ndof = block.Size(seq);

    //std::cout << "posi = "<< pos << std::endl;
    TPZFMatrix<REAL> sol = cmesh->Solution();
    if (pos+1 >= sol.Rows() || ndof < 2) return 0.0;

    return sol(pos+1,0);
}

#include <fstream>
#include <iomanip>
#include <vector> // <-- ADICIONE ISTO
// Eq. (4.123) – passo preditor (k = 1)
static STATE compute_dlambda0_riks(const TPZFMatrix<STATE>& dwb,
                                   const TPZFMatrix<STATE>& dw,
                                   STATE L)
{
    const STATE s    = Dot(dw, dwb);                 // Δu^T * dū
    const STATE ndwb = Norm(dwb);
    const STATE signum = (s > 0.0 ? -1.0 : 1.0);     // Souza Neto 4.123
    return signum * L / ndwb;
}

// Eqs. (4.116) + (4.118) – escolha da raiz para k > 1
static STATE compute_dlambda_riks(const TPZFMatrix<STATE>& dwb,
                                  const TPZFMatrix<STATE>& dws,
                                  const TPZFMatrix<STATE>& dw,
                                  STATE L, int& rootIdx)
{
    const STATE aa = Dot(dwb, dwb);
    TPZFMatrix<STATE> t = dw;
    t += dws;               // t = dw + dws
    const STATE bb = 2.0 * Dot(dwb, t);
    const STATE cc = Dot(t,t) - L*L;

    const STATE eps = 1e-14;
    if (aa < eps) {                                   // cai para linear
        rootIdx = 1;
        return (std::fabs(bb) > eps) ? (-cc/bb) : 0.0;
    }

    STATE disc = bb*bb - 4.0*aa*cc;
    if (disc < 0.0) disc = 0.0;                       // clamp numérico
    const STATE sq = std::sqrt(disc);

    const STATE dl1 = (-bb - sq) / (2.0*aa);          // “menor”
    const STATE dl2 = (-bb + sq) / (2.0*aa);          // “maior”

    auto score = [&](STATE dl)->STATE {
        TPZFMatrix<STATE> x = dw;                     // Δu^(k-1)
        TPZFMatrix<STATE> tmp = dwb; tmp *= dl;
        x += dws; x += tmp;                           // Δu^(k-1)+δu*+δλ dū
        return Dot(x, dw);                            // maximiza (4.118)
    };
    const STATE s1 = score(dl1), s2 = score(dl2);
    if (s1 > s2) { rootIdx = 1; return dl1; }
    else         { rootIdx = 2; return dl2; }
}
REAL IterativeProcessArcLength2(TPZElastoPlasticAnalysis &an,
                                int nsteps,
                                STATE lambda0,
                                STATE L0,
                                std::string vtkfile,int matid,STATE x,STATE y)
{
    // ======== ARQUIVOS DE SAÍDA (reuso dos nomes p/ Python) ========
    std::ofstream out_ld("arc_load_displacement.txt");
    out_ld << std::scientific << std::setprecision(15);
    out_ld << "# step lambda uy\n";

    std::ofstream out_res("arc_residuals.txt");
    out_res << std::scientific << std::setprecision(15);
    out_res << "# step iter lambda L(reserved) res_norm\n";

    auto cmesh = an.Mesh();
    cmesh->Solution().Zero();

    auto* bodymat = dynamic_cast<TMatElastoPlaticMC*>(cmesh->FindMaterial(matid));
    auto* bcmat = dynamic_cast<TPZBndCondT<STATE>*>(cmesh->FindMaterial(matid));

    auto* matmem = dynamic_cast<TPZMatWithMem<TPZElastoPlasticMem>*>(cmesh->FindMaterial(1));
    matmem->ResetMemory();

    TPZManVector<REAL,3> f0;
    if(bodymat)
    {
        std::cout << "material de volume ENCONTRADO \n";
    }else
    {
        if(!bcmat)
        {
            std::cout << "material de contorno nao encontrado \n";
            DebugStop();
        }
        std::cout << "material de contorno ENCONTRADO \n";

    }
    REAL uy=0.;

    TPZFMatrix<STATE> FEXT;

    if(bodymat)
    {
        bodymat->SetBodyForce(bodymat->GetBodyForce0());
        an.Assemble();
        FEXT = an.Rhs();
        TPZManVector<REAL,3> f0={0,0,0};
        bodymat->SetBodyForce(f0);
    }else
    {
        if(!bcmat)
        {
            std::cout << "material de contorno nao encontrado \n";
            DebugStop();
        }
        an.Assemble();
        FEXT = an.Rhs();
        bcmat->Val2()[0]=0.;
        bcmat->Val2()[1]=0.;
        bcmat->Val2()[2]=0.;
    }

    cmesh->Solution().Zero();
    int dim = cmesh->Dimension();
    REAL lambda = lambda0;   // lambda corrente (pode ser alterado durante os passos)
    REAL lambdan = lambda;   // lambda aceito (para rollback seguro)
    REAL L = L0;             // comprimento de arco atual

    // Solução acumulada aceita (u_acc) e incremento acumulado do passo (dw)
    TPZFMatrix<STATE> u_acc = an.Solution();
    u_acc.Zero();

    int step = 0;
    const int max_cuts = 12;
    int cut_count = 0;

    while (step < nsteps) {

        // Estado de trabalho parte do aceito
        TPZFMatrix<STATE> u  = u_acc;
        bool okconv = false;
        REAL normR  = 1e30;
        TPZFMatrix<STATE> dw = an.Solution();
        dw.Zero();

        {
            const int  maxit_inner = 20;
            const REAL etol_inner  = 1e-8;

            TPZFMatrix<STATE> rhs, R, dws, dwb;

            int it = 0;
            int rootIdx = 0;

            // ----- LOG de resíduos deste step (armazenado, só imprime se convergir)
            struct ResEntry { int iter; STATE lambda; STATE Ldummy; STATE res; };
            std::vector<ResEntry> res_log; // <-- NOVO

            std::cout << " \n step = "<< step <<"\n";
            while (it < maxit_inner && normR > etol_inner)
            {
                an.LoadSolution(u);
                an.Assemble();
                rhs = an.Rhs();
                R = FEXT * lambda;
                R += rhs;

                an.Rhs() = R;
                an.Solve();
                dws = an.Solution();

                an.Rhs() = FEXT;
                an.Solve();
                dwb = an.Solution();

                REAL dl = (it == 0) ? compute_dlambda0_riks(dwb, dw, L)
                : compute_dlambda_riks(dwb, dws, dw, L, /*rootIdx*/ *(int[]){0});

                const REAL dl_max = 1;
                if (dl >  dl_max) dl =  dl_max;
                if (dl < -dl_max) dl = -dl_max;

                TPZFMatrix<STATE> dwtot = dwb*dl + dws;

                u      += dwtot;
                dw     += dwtot;
                lambda += dl;

                an.LoadSolution(u);
                an.Assemble();

                rhs   = an.Rhs();
                R     = FEXT * lambda;
                R    += rhs;
                normR = Norm(R);

                STATE normRFEXT = Norm(R)/Norm(FEXT);
                // ---- Guarda o residual desta iteração (só será impresso se o step convergir)
                res_log.push_back({it, lambda, L, normR}); // <-- NOVO

                cout << " iter = " << it
                << "  lambda = " << lambda
                << " dl = "     << dl
                << "  L = "     << L
                //<< " Norm(dw) = "  << Norm(dwtot)
                //<< " normRFEXT = "  << normRFEXT
                << " normR = "  << normR << "\n";
                ++it;

                // if (normR > 1e5) {
                //     break;
                // }
            }

            okconv = (normR <= etol_inner);

            if (okconv) {
                // ---- imprime resíduos deste step (apenas em caso de CONVERGÊNCIA)
                for (const auto &r : res_log) {
                    out_res << step      << " "
                    << (r.iter+0) << " "
                    << r.lambda   << " "
                    << r.Ldummy   << " "
                    << r.res      << "\n";
                }
            }
        }

        if (okconv) {
            // aceita o passo
            cmesh->LoadSolution(u);
            uy+=  UyAtNode(cmesh, x,y);
            out_ld << step   << " " << lambda << " " << uy << "\n";
            an.AcceptSolution();
            u_acc   = an.Solution();

            // salva saída VTK (função externa no seu projeto)
            PostElastoplastic(cmesh, vtkfile, /*matid*/1, /*out_step*/ step , dim);

            // pronto para o próximo step
            ++step;
            cut_count = 0;

            // atualiza “lambda aceito” para rollback seguro
            if(fabs(lambda-lambdan)<1.e-3)break;
            lambdan = lambda;
        } else {
            // rollback e corta L
            an.LoadSolution(u_acc);
            lambda = lambdan;   // volta lambda aceito
            L *= 0.5;
            ++cut_count;

            if (cut_count > 12 || L < 1e-14) {
                std::cout << "[ArcLength] Falha em convergir no step " << step
                << " após " << cut_count << " cortes de L. Abortando ciclo.\n";
                break; // evita loop infinito do ciclo
            }

            // tenta novamente o MESMO step com L menor
            continue;
        }

    } // fim while step < nsteps

    an.AcceptSolution();
    return 0.0;
}
