
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

//#include "TPZPlasticityTest.h"
#include <iostream>
#include <cstdlib>
#include "pzelastoplastic.h"
#include "pzporous.h"
#include "TPZThermoForceA.h"
#include "TPZElasticResponse.h"
#include "pzelastoplasticanalysis.h"
#include "pzanalysis.h"
#include "pzskylstrmatrix.h"
#include "TPZTensor.h"
#include "pzcompelpostproc.h"
#include "pzpostprocmat.h"
#include "pzpostprocanalysis.h"
#include "TPZYCVonMises.h"
//#include "TPZVonMises.h"
#include "pzfstrmatrix.h"
#include "pzbndmat.h"
#include "pzgeoquad.h"
#include "TPZGeoCube.h"
#include "pzgeotetrahedra.h"
#include "pzgeopyramid.h"
#include "tpzgeoelrefpattern.h"
#include "pzbndcond.h"
#include "pzstepsolver.h"
#include "TPZTensor.h"
#include "TPZYCMohrCoulomb.h"
#include "TPZMohrCoulomb.h"
//#include "TPZVonMises.h"
//#include "TPZYCVonMises.h"
#include "TPZDruckerPrager.h"
#include "GeoMeshClass.h"
#include "pzelastoplastic2D.h"
#include <pzmathyperelastic.h>
#include "tpzycvonmisescombtresca.h"
#include "TPZMohrCoulombNeto.h"
#include "TPZSandlerDimaggio.h"
#include "clock_timer.h"
#include "TPBrAcidFunc.h"


#include "pzcmesh.h"
#include "TPZSandlerDimaggio.h"
#include "TPZTensor.h"
#include "pzgeoel.h"
#include "pzpostprocanalysis.h"
#include "pzsandlerextPV.h"
#include "TPZPlasticStepPV.h"
#include "pzstack.h"
#include "TPZYCMohrCoulombPV.h"
#include "pzelasticSest2D.h"
#include "pzelastoplasticSest2D.h"
#include "TPBrAcidFunc.h"
#include "TPZElasticCriteria.h"
#include "pzstring.h"


using namespace pzshape; // needed for TPZShapeCube and related classes



#include "pzlog.h"
//#include "tpztimer.h"
#include "TPZTimer.h"
#include "WellBoreAnalysis.h"
#include "pzbfilestream.h"
#include "TPZProjectEllipse.h"
#include "arglib.h"
#include "run_stats_table.h"

#define MACOS
#ifdef MACOS

#include <iostream>
#include <math.h>
#include <signal.h>
#include <fenv.h>
#include <xmmintrin.h>

#define ENABLE_FPO_EXCEPTIONS _MM_SET_EXCEPTION_MASK(_MM_GET_EXCEPTION_MASK() & ~_MM_MASK_INVALID);


#define DECLARE_FPO_HANDLER_FUNC void InvalidFPOHandler(int signo) {\
switch(signo) {\
case SIGFPE: std::cout << "ERROR : Invalid Arithmetic operation." << std::endl; break;\
}\
exit(signo);\
}

#define ATTACH_FPO_SIGNAL struct sigaction act = {};\
act.sa_handler = InvalidFPOHandler;\
sigaction(SIGFPE, &act, NULL);


DECLARE_FPO_HANDLER_FUNC;
#endif

#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("plasticity.main"));
#endif

#ifdef LOG4CXX
static LoggerPtr loggerEllipse(Logger::getLogger("LogEllipse"));
#endif

#ifdef USING_TBB
#include "tbb/task_scheduler_init.h"
using namespace tbb;
// If you have issues with: dyld: Library not loaded: libtbb.dylib
// try setting the LD path. Ex:
//   export DYLD_FALLBACK_LIBRARY_PATH=/Users/borin/Desktop/neopz/tbb40_297oss/lib/
#endif


RunStatsTable plast_tot("-tpz_plast_tot", "Raw data table statistics for the main execution.");
clarg::argInt NumberOfThreads("-nt", "Number of threads for WellBoreAnalysis", 1);

TPZGeoMesh * ShellGeoGen(int ndiv, REAL alpha, REAL Rm, REAL h, REAL a, REAL b, int ndirectdiv);
void RefinamentoUniforme(TPZGeoMesh & gMesh, int &nh);
void Config1()
{
    //EVertical
    //ENonPenetrating
    InitializePZLOG();
    gRefDBase.InitializeUniformRefPattern(EOned);
    gRefDBase.InitializeUniformRefPattern(ETriangle);
    gRefDBase.InitializeUniformRefPattern(EQuadrilateral);
    std::cout << std::setprecision(15);
    
    TPZWellBoreAnalysis well;
    
#ifndef USING_TBB
    if (NumberOfThreads.get_value()>=0) {
        TPZWellBoreAnalysis::TConfig::gNumThreads=NumberOfThreads.get_value();
    }
#else
    int number_tbb=NumberOfThreads.get_value();
    if(number_tbb<=0)number_tbb=1;
    task_scheduler_init init(number_tbb);
#endif

    
    STATE biotcoef = 0.659;
    well.SetBiotCoefficient(biotcoef);
    
    REAL reservoirPressure=57.2;
    well.SetReservoirPressure(reservoirPressure);
    
    
    REAL SH,Sh,SV;
    Sh=-83.5;
    SH=-99.8;
    SV=-85.9;
    TPZManVector<STATE,3> confinementTotal(3,0.);
    confinementTotal[0] = Sh;
    confinementTotal[1] = SH;
    confinementTotal[2] = SV;
    REAL WellPressure = 57.2; //66.6 61.1 57.2
    well.SetConfinementTotalStresses(confinementTotal, WellPressure);
    
    
    REAL innerradius = 4.25*0.0254;
    REAL outerradius = 3.;
    well.SetInnerOuterRadius(innerradius, outerradius);
    
    REAL area,removed;
    std::string output = "Config1.vtk";
    well.SetVtkOutPutName(output);
 

    REAL sqj2_refine=1.e-9;
    int Startfrom=0;
    const int nsubsteps = 5;
    if (Startfrom == 0)
    {
        well.SetInnerOuterRadius(innerradius, outerradius);
        
        REAL poisson = 0.203;
        REAL elast = 29269.;
        REAL A = 152.54;
        REAL B = 0.0015489;
        REAL C = 146.29;
        REAL R = 0.91969;
        REAL D = 0.018768;
        REAL W = 0.006605;
        
        bool modelMC =false;
        
        if (modelMC)
        {
            REAL cohesion = 13.;
            REAL Phi = 0.52;
            well.SetMohrCoulombParameters(poisson, elast, cohesion, Phi, Phi);
            
        }
        else
        {
            well.SetSanderDiMaggioParameters(poisson, elast, A, B, C, R, D, W);
            
            
        }
        
        int porder = 2;
        int nrad=20;
        int ncircle = 40;
        REAL delx = 0.5*innerradius*M_PI_2/ncircle;
        TPZManVector<int,2> numdiv(2);
        numdiv[0] = nrad;
        numdiv[1] = ncircle;
        well.SetMeshTopology(delx, numdiv);
        
        
        well.GetCurrentConfig()->fWellConfig = EVerticalWell;
        
        
        well.GetCurrentConfig()->CreateGeometricMesh();
        
        
        well.GetCurrentConfig()->CreateComputationalMesh(porder);
        
        
        well.GetCurrentConfig()->CreatePostProcessingMesh();

        well.PostProcess(0);
        
    }
    if (Startfrom ==0)
    {
        
        int nsteps = 5;
        int numnewton = 90;
        well.ExecuteInitialSimulation(nsteps, numnewton);

        area = well.GetCurrentConfig()->ComputeAreaAboveSqJ2(1.e-9);
        removed = well.GetCurrentConfig()->RemovedArea();
        
        cout  << " plastified area initial = "<< area << endl;
        cout  << " removed area initial = "<< removed << endl;
        
        well.PostProcess(0);
        
        TPZBFileStream save;
        save.OpenWrite("Config1-0.bin");
        well.Write(save);
        
    }
    
    if (Startfrom ==1)
    {
        TPZBFileStream read;
        read.OpenRead("Config1-0.bin");
        well.Read(read);
    }
    
    
    if (Startfrom <=1)
    {
        std::cout << "\n ------- 1 -------- "<<std::endl;
        
        std::stringstream out;
        out << "First pass\n";
        well.PRefineElementAbove(sqj2_refine, 3,out);
        well.DivideElementsAbove(sqj2_refine,out);
        
        {
            std::ofstream out("gmesh_good.txt");
            well.GetCurrentConfig()->fGMesh.Print(out);
            std::ofstream co("cmesh_good.txt");
            well.GetCurrentConfig()->fCMesh.Print(co);
        }
        
        
        int substepsloc  = 1;
        well.ExecuteSimulation(substepsloc,out);
        
        area = well.GetCurrentConfig()->ComputeAreaAboveSqJ2(1.e-9);
        removed = well.GetCurrentConfig()->RemovedArea();
        
        cout  << " plastified area step 1 = "<< area << endl;
        cout  << " removed area step 1 = "<< removed << endl;
        
        well.PostProcess(0);
        REAL a, b;
        well.ComputeAandB(sqj2_refine, a,b);
        if(a >innerradius )
        {
            well.AddEllipticBreakout(a, b,out);
            int substepsloc = 1;
            well.ExecuteSimulation(substepsloc,out);
            
            area = well.GetCurrentConfig()->ComputeAreaAboveSqJ2(1.e-9);
            removed = well.GetCurrentConfig()->RemovedArea();
            
            
            cout  << " plastified area  step 1 after addElliptc = "<< area << endl;
            cout  << " removed area  step 1 after addElliptc = "<< removed << endl;
            
            
            well.PostProcess(0);
        }
        TPZBFileStream save;
        save.OpenWrite("Config1-1.bin");
        well.Write(save);
        
        
    }
    
    if (Startfrom ==2)
    {
        TPZBFileStream read;
        read.OpenRead("Config1-1.bin");
        well.Read(read);
    }
    
    if (Startfrom <=2)
    {
        
        std::cout << "\n ------- 2 -------- "<<std::endl;
        std::stringstream out;
        out << "Applying second pass\n";
        well.PRefineElementAbove(sqj2_refine, 3,out);
        well.DivideElementsAbove(sqj2_refine,out);
        int substepsloc = 1;
        well.ExecuteSimulation(substepsloc,out);
        
        
        area = well.GetCurrentConfig()->ComputeAreaAboveSqJ2(1.e-9);
        removed = well.GetCurrentConfig()->RemovedArea();
        cout  << " plastified area step 2 = "<< area << endl;
        cout  << " removed area step 2 = "<< removed << endl;
        
        
        well.PostProcess(0);
        REAL a, b;
        well.ComputeAandB(sqj2_refine, a,b);
        if(a >innerradius )
        {
            well.AddEllipticBreakout(a, b,out);
            int substepsloc = 1;
            well.ExecuteSimulation(substepsloc,out);
            
            
            area = well.GetCurrentConfig()->ComputeAreaAboveSqJ2(1.e-9);
            removed = well.GetCurrentConfig()->RemovedArea();
            
            cout  << " plastified area  step 2 after addElliptc  = "<< area << endl;
            cout  << " removed area  step 2 after addElliptc = "<< removed << endl;
            
            well.PostProcess(0);
        }
        well.AppendExecutionLog(out);
        TPZBFileStream save;
        save.OpenWrite("Config1-2.bin");
        well.Write(save);
        
    }
    
    
    
    if (Startfrom ==3)
    {
        TPZBFileStream read;
        read.OpenRead("Config1-2.bin");
        well.Read(read);
    }
    
    if (Startfrom <=3)
    {
        
        std::cout << "\n ------- 3 -------- "<<std::endl;
        std::stringstream out;
        out << "Applying third pass\n";
        well.PRefineElementAbove(sqj2_refine, 3,out);
        well.DivideElementsAbove(sqj2_refine,out);
        int substepsloc =1;
        well.ExecuteSimulation(substepsloc,out);
        
        area = well.GetCurrentConfig()->ComputeAreaAboveSqJ2(1.e-9);
        removed = well.GetCurrentConfig()->RemovedArea();
        cout  << " plastified area step 3 = "<< area << endl;
        cout  << " removed area step 3 = "<< removed << endl;
        
        well.PostProcess(0);
        REAL a, b;
        well.ComputeAandB(sqj2_refine, a,b);
        if(a >innerradius )
        {
            well.AddEllipticBreakout(a, b,out);
            well.ExecuteSimulation(substepsloc,out);
            
            
            area = well.GetCurrentConfig()->ComputeAreaAboveSqJ2(1.e-9);
            removed = well.GetCurrentConfig()->RemovedArea();
            
            cout  << " plastified area  step 3 after addElliptc  = "<< area << endl;
            cout  << " removed area  step 3 after addElliptc = "<< removed << endl;
            
            well.PostProcess(0);
        }
        well.AppendExecutionLog(out);
        TPZBFileStream save;
        save.OpenWrite("Config1-3.bin");
        well.Write(save);
        
        
    }
    
    
    
    if (Startfrom ==4)
    {
        TPZBFileStream read;
        read.OpenRead("Config1-3.bin");
        well.Read(read);
    }
    
    if (Startfrom <=4)
    {
        
        std::cout << "\n ------- 4 -------- "<<std::endl;
        std::stringstream out;
        out << "Applying fourth pass\n";
        well.PRefineElementAbove(sqj2_refine, 3,out);
        well.DivideElementsAbove(sqj2_refine,out);
        int substepsloc =1;
        well.ExecuteSimulation(substepsloc,out);
        well.PostProcess(0);
        REAL a, b;
        well.ComputeAandB(sqj2_refine, a,b);
        if(a >innerradius )
        {
            well.AddEllipticBreakout(a, b,out);
            well.ExecuteSimulation(substepsloc,out);
            well.PostProcess(0);
        }
        well.AppendExecutionLog(out);
        TPZBFileStream save;
        save.OpenWrite("Config1-4.bin");
        well.Write(save);
        
        
    }
    well.PrintExecutionLog(std::cout);
}


void Config2()
{
    //EHorizontalWellalongH
    //ENonPenetrating
    InitializePZLOG();
    gRefDBase.InitializeUniformRefPattern(EOned);
    gRefDBase.InitializeUniformRefPattern(ETriangle);
    gRefDBase.InitializeUniformRefPattern(EQuadrilateral);
    std::cout << std::setprecision(15);
    
    TPZWellBoreAnalysis well;
    
    STATE biotcoef = 0.659;
    well.SetBiotCoefficient(biotcoef);
    
    REAL reservoirPressure=57.2;
    well.SetReservoirPressure(reservoirPressure);
    
    
    REAL SH,Sh,SV;
    Sh=-83.5;
    SH=-99.8;
    SV=-85.9;
    TPZManVector<STATE,3> confinementTotal(3,0.);
    confinementTotal[0] = Sh;
    confinementTotal[1] = SH;
    confinementTotal[2] = SV;
    REAL WellPressure = 57.2; //66.6 61.1 57.2
    well.SetConfinementTotalStresses(confinementTotal, WellPressure);
    
    
    REAL innerradius = 4.25*0.0254;
    REAL outerradius = 3.;
    well.SetInnerOuterRadius(innerradius, outerradius);
    
    
    std::string output = "Config2.vtk";
    well.SetVtkOutPutName(output);
    
    
    REAL sqj2_refine=0.0001;
    const int nsubsteps = 5;
    int Startfrom=0;
    if (Startfrom == 0)
    {
        well.SetInnerOuterRadius(innerradius, outerradius);
        
        REAL poisson = 0.203;
        REAL elast = 29269.;
        REAL A = 152.54;
        REAL B = 0.0015489;
        REAL C = 146.29;
        REAL R = 0.91969;
        REAL D = 0.018768;
        REAL W = 0.006605;
        
        bool modelMC =false;
        
        if (modelMC)
        {
            REAL cohesion = 13.;
            REAL Phi = 0.52;
            well.SetMohrCoulombParameters(poisson, elast, cohesion, Phi, Phi);
            
        }
        else
        {
            well.SetSanderDiMaggioParameters(poisson, elast, A, B, C, R, D, W);
            
            
        }
        
        int porder = 2;
        int nrad=20;
        int ncircle = 40;
        REAL delx = 0.5*innerradius*M_PI_2/ncircle;
        TPZManVector<int,2> numdiv(2);
        numdiv[0] = nrad;
        numdiv[1] = ncircle;
        well.SetMeshTopology(delx, numdiv);
        
        
        well.GetCurrentConfig()->fWellConfig = EHorizontalWellalongH;
        
        
        well.GetCurrentConfig()->CreateGeometricMesh();
        
        
        well.GetCurrentConfig()->CreateComputationalMesh(porder);
        
        
        well.GetCurrentConfig()->CreatePostProcessingMesh();
        
        well.PostProcess(0);
        
    }
    if (Startfrom ==0)
    {
        
        int nsteps = 5;
        int numnewton = 90;
        well.GetCurrentConfig()->ModifyWellElementsToQuadratic();
        well.ExecuteInitialSimulation(nsteps, numnewton);
        well.PostProcess(0);
        TPZBFileStream save;
        save.OpenWrite("Config2-0.bin");
        well.Write(save);
        
    }
    
    if (Startfrom ==1)
    {
        TPZBFileStream read;
        read.OpenRead("Config2-0.bin");
        well.Read(read);
    }
    
    
    if (Startfrom <=1)
    {
        std::cout << "\n ------- 1 -------- "<<std::endl;
        std::stringstream out;
        out << "Applying 1st pass\n";

        well.PRefineElementAbove(sqj2_refine, 3,out);
        well.DivideElementsAbove(sqj2_refine,out);
        well.ExecuteSimulation(nsubsteps,out);
        well.PostProcess(0);
        REAL a, b;
        well.ComputeAandB(sqj2_refine, a,b);
        if(a >innerradius )
        {
            well.AddEllipticBreakout(a, b,out);
            well.ExecuteSimulation(nsubsteps,out);
            well.PostProcess(0);
        }
        TPZBFileStream save;
        save.OpenWrite("Config2-1.bin");
        well.Write(save);
        
        
    }
    
    if (Startfrom ==2)
    {
        TPZBFileStream read;
        read.OpenRead("Config2-1.bin");
        well.Read(read);
    }
    
    if (Startfrom <=2)
    {
        std::stringstream out;
        out << "Applying 2nd pass\n";
        
        std::cout << "\n ------- 2 -------- "<<std::endl;
        well.PRefineElementAbove(sqj2_refine, 3,out);
        well.DivideElementsAbove(sqj2_refine,out);
        well.ExecuteSimulation(nsubsteps,out);
        well.PostProcess(0);
        REAL a, b;
        well.ComputeAandB(sqj2_refine, a,b);
        if(a >innerradius )
        {
            well.AddEllipticBreakout(a, b,out);
            well.ExecuteSimulation(nsubsteps,out);
            well.PostProcess(0);
        }
        TPZBFileStream save;
        save.OpenWrite("Config2-2.bin");
        well.Write(save);
        
    }
    
    
    
    if (Startfrom ==3)
    {
        TPZBFileStream read;
        read.OpenRead("Config2-2.bin");
        well.Read(read);
    }
    
    if (Startfrom <=3)
    {
        
        std::cout << "\n ------- 3 -------- "<<std::endl;
        std::stringstream out;
        out << "Applying 3nd pass\n";
        well.PRefineElementAbove(sqj2_refine, 3,out);
        well.DivideElementsAbove(sqj2_refine,out);
        well.ExecuteSimulation(nsubsteps,out);
        well.PostProcess(0);
        REAL a, b;
        well.ComputeAandB(sqj2_refine, a,b);
        if(a >innerradius )
        {
            well.AddEllipticBreakout(a, b,out);
            well.ExecuteSimulation(nsubsteps,out);
            well.PostProcess(0);
        }
        TPZBFileStream save;
        save.OpenWrite("Config2-3.bin");
        well.Write(save);
        
        
    }
    
    
    
    if (Startfrom ==4)
    {
        TPZBFileStream read;
        read.OpenRead("Config2-3.bin");
        well.Read(read);
    }
    
    if (Startfrom <=4)
    {
        
        std::cout << "\n ------- 4 -------- "<<std::endl;
        std::stringstream out;
        out << "Applying 4th pass\n";
        well.PRefineElementAbove(sqj2_refine, 3,out);
        well.DivideElementsAbove(sqj2_refine,out);
        well.ExecuteSimulation(nsubsteps,out);
        well.PostProcess(0);
        REAL a, b;
        well.ComputeAandB(sqj2_refine, a,b);
        if(a >innerradius )
        {
            well.AddEllipticBreakout(a, b,out);
            well.ExecuteSimulation(nsubsteps,out);
            well.PostProcess(0);
        }
        TPZBFileStream save;
        save.OpenWrite("Config2-4.bin");
        well.Write(save);
        
        
    }
    
    if (Startfrom ==5)
    {
        TPZBFileStream read;
        read.OpenRead("Config2-4.bin");
        well.Read(read);
    }
    
    if (Startfrom <=5)
    {
        
        std::cout << "\n ------- 5 -------- "<<std::endl;
        std::stringstream out;
        out << "Applying 5th pass\n";
        well.PRefineElementAbove(sqj2_refine, 3,out);
        well.DivideElementsAbove(sqj2_refine,out);
        well.ExecuteSimulation(nsubsteps,out);
        well.PostProcess(0);
        REAL a, b;
        well.ComputeAandB(sqj2_refine, a,b);
        if(a >innerradius )
        {
            well.AddEllipticBreakout(a, b,out);
            well.ExecuteSimulation(nsubsteps,out);
            well.PostProcess(0);
        }
        TPZBFileStream save;
        save.OpenWrite("Config2-5.bin");
        well.Write(save);
        
        
    }
    
}


void Config3()
{
    //EHorizontalWellalongh
    //ENonPenetrating
    InitializePZLOG();
    gRefDBase.InitializeUniformRefPattern(EOned);
    gRefDBase.InitializeUniformRefPattern(ETriangle);
    gRefDBase.InitializeUniformRefPattern(EQuadrilateral);
    std::cout << std::setprecision(15);
    
    TPZWellBoreAnalysis well;
    
    STATE biotcoef = 0.659;
    well.SetBiotCoefficient(biotcoef);
    
    REAL reservoirPressure=57.2;
    well.SetReservoirPressure(reservoirPressure);
    
    
    REAL SH,Sh,SV;
    Sh=-83.5;
    SH=-99.8;
    SV=-85.9;
    TPZManVector<STATE,3> confinementTotal(3,0.);
    confinementTotal[0] = Sh;
    confinementTotal[1] = SH;
    confinementTotal[2] = SV;
    REAL WellPressure = 57.2; //66.6 61.1 57.2
    well.SetConfinementTotalStresses(confinementTotal, WellPressure);
    
    
    REAL innerradius = 4.25*0.0254;
    REAL outerradius = 3.;
    well.SetInnerOuterRadius(innerradius, outerradius);
    
    
    std::string output = "Config3.vtk";
    well.SetVtkOutPutName(output);
    
    
    REAL sqj2_refine=0.0001;
    const int nsubsteps = 5;
    int Startfrom=0;
    if (Startfrom == 0)
    {
        well.SetInnerOuterRadius(innerradius, outerradius);
        
        REAL poisson = 0.203;
        REAL elast = 29269.;
        REAL A = 152.54;
        REAL B = 0.0015489;
        REAL C = 146.29;
        REAL R = 0.91969;
        REAL D = 0.018768;
        REAL W = 0.006605;
        
        bool modelMC =false;
        
        if (modelMC)
        {
            REAL cohesion = 13.;
            REAL Phi = 0.52;
            well.SetMohrCoulombParameters(poisson, elast, cohesion, Phi, Phi);
            
        }
        else
        {
            well.SetSanderDiMaggioParameters(poisson, elast, A, B, C, R, D, W);
            
            
        }
        
        int porder = 2;
        int nrad=20;
        int ncircle = 40;
        REAL delx = 0.5*innerradius*M_PI_2/ncircle;
        TPZManVector<int,2> numdiv(2);
        numdiv[0] = nrad;
        numdiv[1] = ncircle;
        well.SetMeshTopology(delx, numdiv);
        
        
        well.GetCurrentConfig()->fWellConfig = EHorizontalWellalongh;
        
        
        well.GetCurrentConfig()->CreateGeometricMesh();
        
        
        well.GetCurrentConfig()->CreateComputationalMesh(porder);
        
        
        well.GetCurrentConfig()->CreatePostProcessingMesh();
        
        well.PostProcess(0);
        
    }
    if (Startfrom ==0)
    {
        
        int nsteps = 5;
        int numnewton = 90;
        well.GetCurrentConfig()->ModifyWellElementsToQuadratic();
        well.ExecuteInitialSimulation(nsteps, numnewton);
        well.PostProcess(0);
        TPZBFileStream save;
        save.OpenWrite("Config3-0.bin");
        well.Write(save);
        
    }
    
    if (Startfrom ==1)
    {
        TPZBFileStream read;
        read.OpenRead("Config3-0.bin");
        well.Read(read);
    }
    
    
    if (Startfrom <=1)
    {
        std::cout << "\n ------- 1 -------- "<<std::endl;
        
        std::stringstream out;
        out << "Applying 1st pass\n";
        well.PRefineElementAbove(sqj2_refine, 3,out);
        well.DivideElementsAbove(sqj2_refine,out);
        well.ExecuteSimulation(nsubsteps,out);
        well.PostProcess(0);
        REAL a, b;
        well.ComputeAandB(sqj2_refine, a,b);
        if(a >innerradius )
        {
            well.AddEllipticBreakout(a, b,out);
            well.ExecuteSimulation(nsubsteps,out);
            well.PostProcess(0);
        }
        TPZBFileStream save;
        save.OpenWrite("Config3-1.bin");
        well.Write(save);
        
        
    }
    
    if (Startfrom ==2)
    {
        TPZBFileStream read;
        read.OpenRead("Config3-1.bin");
        well.Read(read);
    }
    
    if (Startfrom <=2)
    {
        std::stringstream out;
        out << "Applying 2nd pass\n";
        
        std::cout << "\n ------- 2 -------- "<<std::endl;
        well.PRefineElementAbove(sqj2_refine, 3,out);
        well.DivideElementsAbove(sqj2_refine,out);
        well.ExecuteSimulation(nsubsteps,out);
        well.PostProcess(0);
        REAL a, b;
        well.ComputeAandB(sqj2_refine, a,b);
        if(a >innerradius )
        {
            well.AddEllipticBreakout(a, b,out);
            well.ExecuteSimulation(nsubsteps,out);
            well.PostProcess(0);
        }
        TPZBFileStream save;
        save.OpenWrite("Config3-2.bin");
        well.Write(save);
        
    }
    
    
    
    if (Startfrom ==3)
    {
        TPZBFileStream read;
        read.OpenRead("Config3-2.bin");
        well.Read(read);
    }
    
    if (Startfrom <=3)
    {
        
        std::cout << "\n ------- 3 -------- "<<std::endl;
        std::stringstream out;
        out << "Applying 3rd pass\n";
        well.PRefineElementAbove(sqj2_refine, 3,out);
        well.DivideElementsAbove(sqj2_refine,out);
        well.ExecuteSimulation(nsubsteps,out);
        well.PostProcess(0);
        REAL a, b;
        well.ComputeAandB(sqj2_refine, a,b);
        if(a >innerradius )
        {
            well.AddEllipticBreakout(a, b,out);
            well.ExecuteSimulation(nsubsteps,out);
            well.PostProcess(0);
        }
        TPZBFileStream save;
        save.OpenWrite("Config3-3.bin");
        well.Write(save);
        
        
    }
    
    
    
    if (Startfrom ==4)
    {
        TPZBFileStream read;
        read.OpenRead("Config3-3.bin");
        well.Read(read);
    }
    
    if (Startfrom <=4)
    {
        std::stringstream out;
        out << "Applying 4th pass\n";
        
        std::cout << "\n ------- 4 -------- "<<std::endl;
        well.PRefineElementAbove(sqj2_refine, 3,out);
        well.DivideElementsAbove(sqj2_refine,out);
        well.ExecuteSimulation(nsubsteps,out);
        well.PostProcess(0);
        REAL a, b;
        well.ComputeAandB(sqj2_refine, a,b);
        if(a >innerradius )
        {
            well.AddEllipticBreakout(a, b,out);
            well.ExecuteSimulation(nsubsteps,out);
            well.PostProcess(0);
        }
        TPZBFileStream save;
        save.OpenWrite("Config3-4.bin");
        well.Write(save);
        
        
    }
    
    if (Startfrom ==5)
    {
        TPZBFileStream read;
        read.OpenRead("Config3-4.bin");
        well.Read(read);
    }
    
    if (Startfrom <=5)
    {
        
        std::stringstream out;
        out << "Applying 5th pass\n";
        std::cout << "\n ------- 5 -------- "<<std::endl;
        well.PRefineElementAbove(sqj2_refine, 3,out);
        well.DivideElementsAbove(sqj2_refine,out);
        well.ExecuteSimulation(nsubsteps,out);
        well.PostProcess(0);
        REAL a, b;
        well.ComputeAandB(sqj2_refine, a,b);
        if(a >innerradius )
        {
            well.AddEllipticBreakout(a, b,out);
            well.ExecuteSimulation(nsubsteps,out);
            well.PostProcess(0);
        }

        TPZBFileStream save;
        save.OpenWrite("Config3-5.bin");
        well.Write(save);
        
        
    }
    
}

void Config4()
{
    InitializePZLOG();
    gRefDBase.InitializeUniformRefPattern(EOned);
    gRefDBase.InitializeUniformRefPattern(ETriangle);
    gRefDBase.InitializeUniformRefPattern(EQuadrilateral);
    std::cout << std::setprecision(15);
    
    TPZWellBoreAnalysis well;
    
    STATE biotcoef = 0.659;
    well.SetBiotCoefficient(biotcoef);
    
    REAL reservoirPressure=57.2;
    well.SetReservoirPressure(reservoirPressure);
    
    
    REAL SH,Sh,SV;
    Sh=-83.5;
    SH=-99.8;
    SV=-85.9;
    TPZManVector<STATE,3> confinementTotal(3,0.);
    confinementTotal[0] = Sh;
    confinementTotal[1] = SH;
    confinementTotal[2] = SV;
    REAL WellPressure = 57.2; //66.6 61.1 57.2
    well.SetConfinementTotalStresses(confinementTotal, WellPressure);
    
    
    REAL innerradius = 4.25*0.0254;
    REAL outerradius = 3.;
    well.SetInnerOuterRadius(innerradius, outerradius);
    
    
    std::string output = "Config4.vtk";
    const int nsubsteps = 5;
    well.SetVtkOutPutName(output);
    EPlasticModel Emodel = EElastic;
    if (Emodel == EMohrCoulomb)
    {
        REAL poisson = 0.203;
        REAL elast = 29269.;
        REAL cohesion = 13.;
        REAL Phi = 0.52;
        well.SetMohrCoulombParameters(poisson, elast, cohesion, Phi, Phi);        
    }
    else if (Emodel == ESandler)
    {
        REAL poisson = 0.203;
        REAL elast = 29269.;
        REAL A = 152.54;
        REAL B = 0.0015489;
        REAL C = 146.29;
        REAL R = 0.91969;
        REAL D = 0.018768;
        REAL W = 0.006605;
        well.SetSanderDiMaggioParameters(poisson, elast, A, B, C, R, D, W);
    }
    else if (Emodel == EElastic){ // Mohr-Coulomb with a VERY far way surface
        REAL poisson = 0.203;
        REAL elast = 29269.;
        REAL cohesion = 1.e8; // Very very big
        REAL Phi = 1.5533430342749532; // 89 degrees
        
#ifdef PlasticPQP
      well.SetMohrCoulombParameters(poisson, elast, cohesion, Phi, Phi);
      well.GetCurrentConfig()->fModel = EElastic;
#else
      well.SetElasticParameters(poisson, elast);
      well.GetCurrentConfig()->fModel = EElastic;
#endif
    }
    
    
    int Startfrom=0;
    if (Startfrom == 0)
    {
        int porder = 2;
        int nrad=20;
        int ncircle = 40;
        REAL delx = 0.5*innerradius*M_PI_2/ncircle;
        TPZManVector<int,2> numdiv(2);
        numdiv[0] = nrad;
        numdiv[1] = ncircle;
        well.SetMeshTopology(delx, numdiv);
        well.GetCurrentConfig()->fWellConfig = EVerticalWell;
        well.GetCurrentConfig()->CreateGeometricMesh();
        
        well.GetCurrentConfig()->CreateComputationalMesh(porder);
        std::cout << "Average vertical stress " << well.GetCurrentConfig()->AverageVerticalStress() << std::endl;
        well.GetCurrentConfig()->CreatePostProcessingMesh();
        //REAL farfieldwork = well.GetCurrentConfig()->ComputeFarFieldWork();
        //well.PostProcess(0);
        
    }
    if (Startfrom ==0)
    {
        
        int nsteps = 5;
        int numnewton = 80;
        well.GetCurrentConfig()->ModifyWellElementsToQuadratic();
        well.ExecuteInitialSimulation(nsteps, numnewton);
        well.PostProcess(0);
        TPZBFileStream save;
        save.OpenWrite("wellbore0.bin");
        well.Write(save);
        
    }
    
    if (Startfrom ==1)
    {
        TPZBFileStream read;
        read.OpenRead("wellbore0.bin");
        well.Read(read);
    }
    
    
    if (Startfrom <=1)
    {
        std::stringstream out;
        out << "Applying 1st pass\n";
        std::cout << "\n ------- 1 -------- "<<std::endl;
        
        well.PRefineElementAbove(0.000001, 2,out);
        well.ExecuteSimulation(nsubsteps,out);
        well.PostProcess(0);
        
        cout << "Average vertical stress " << well.GetCurrentConfig()->AverageVerticalStress() << std::endl;
        
        well.SetFluidModel(EPenetrating);
        well.ExecuteSimulation(nsubsteps,out);
        well.PostProcess(0);
        cout << "Penetrating fluid Average vertical stress " << well.GetCurrentConfig()->AverageVerticalStress() << std::endl;
        
        
        well.EvolveBothPressures(4, WellPressure*1.2,reservoirPressure);
        cout << "Higher well pressure Average vertical stress " << well.GetCurrentConfig()->AverageVerticalStress() << std::endl;
        
        well.EvolveBothPressures(4,WellPressure,reservoirPressure*0.8);
        cout << "Lower reservoir pressure Average vertical stress " << well.GetCurrentConfig()->AverageVerticalStress() << std::endl;
        
        TPZBFileStream save;
        save.OpenWrite("wellbore1.bin");
        well.Write(save);
    }
    
}

void Config5()
{
    InitializePZLOG();
    gRefDBase.InitializeUniformRefPattern(EOned);
    gRefDBase.InitializeUniformRefPattern(ETriangle);
    gRefDBase.InitializeUniformRefPattern(EQuadrilateral);
    std::cout << std::setprecision(15);
    
    TPZWellBoreAnalysis well;
    
    STATE biotcoef = 0.659;
    well.SetBiotCoefficient(biotcoef);
    
    REAL reservoirPressure=57.2;
    well.SetReservoirPressure(reservoirPressure);
    
    
    REAL SH,Sh,SV;
    Sh=-83.5;
    SH=-99.8;
    SV=-85.9;
    TPZManVector<STATE,3> confinementTotal(3,0.);
    confinementTotal[0] = Sh;
    confinementTotal[1] = SH;
    confinementTotal[2] = SV;
    REAL WellPressure = 57.2; //66.6 61.1 57.2
    well.SetConfinementTotalStresses(confinementTotal, WellPressure);
    
    
    REAL innerradius = 4.25*0.0254;
    REAL outerradius = 3.;
    well.SetInnerOuterRadius(innerradius, outerradius);
    
    
    std::string output = "Config5.vtk";
    const int nsubsteps = 5;
    well.SetVtkOutPutName(output);
    EPlasticModel Emodel = ESandler;
    if (Emodel == EMohrCoulomb)
    {
        REAL poisson = 0.203;
        REAL elast = 29269.;
        REAL cohesion = 13.;
        REAL Phi = 0.52;
        well.SetMohrCoulombParameters(poisson, elast, cohesion, Phi, Phi);        
    }
    else if (Emodel == ESandler)
    {
        REAL poisson = 0.203;
        REAL elast = 29269.;
        REAL A = 152.54;
        REAL B = 0.0015489;
        REAL C = 146.29;
        REAL R = 0.91969;
        REAL D = 0.018768;
        REAL W = 0.006605;
        well.SetSanderDiMaggioParameters(poisson, elast, A, B, C, R, D, W);
    }
    else if (Emodel == EElastic){ // Mohr-Coulomb with a VERY far way surface
        REAL poisson = 0.203;
        REAL elast = 29269.;
        REAL cohesion = 1.e8; // Very very big
        REAL Phi = 1.5533430342749532; // 89 degrees
        
#ifdef PlasticPQP
      well.SetMohrCoulombParameters(poisson, elast, cohesion, Phi, Phi);
      well.GetCurrentConfig()->fModel = EElastic;
#else
      well.SetElasticParameters(poisson, elast);
      well.GetCurrentConfig()->fModel = EElastic;
#endif
    }
    
    
    int Startfrom=0;
    if (Startfrom == 0)
    {
        int porder = 2;
        int nrad=20;
        int ncircle = 40;
        REAL delx = 0.5*innerradius*M_PI_2/ncircle;
        TPZManVector<int,2> numdiv(2);
        numdiv[0] = nrad;
        numdiv[1] = ncircle;
        well.SetMeshTopology(delx, numdiv);
        well.GetCurrentConfig()->fWellConfig = EVerticalWell;
        well.GetCurrentConfig()->CreateGeometricMesh();
        
        well.GetCurrentConfig()->CreateComputationalMesh(porder);
        std::cout << "Average vertical stress " << well.GetCurrentConfig()->AverageVerticalStress() << std::endl;
        well.GetCurrentConfig()->CreatePostProcessingMesh();
        //REAL farfieldwork = well.GetCurrentConfig()->ComputeFarFieldWork();
        //well.PostProcess(0);
        
    }
    if (Startfrom ==0)
    {
        
        int nsteps = 5;
        int numnewton = 80;
        well.GetCurrentConfig()->ModifyWellElementsToQuadratic();
        well.ExecuteInitialSimulation(nsteps, numnewton);
        well.PostProcess(0);
        TPZBFileStream save;
        save.OpenWrite("wellbore0.bin");
        well.Write(save);
        
    }
    
    if (Startfrom ==1)
    {
        TPZBFileStream read;
        read.OpenRead("wellbore0.bin");
        well.Read(read);
    }
    
    
    if (Startfrom <=1)
    {
        std::stringstream out;
        out << "Applying 1st pass\n";
        std::cout << "\n ------- 1 -------- "<<std::endl;
        
        well.PRefineElementAbove(0.0001, 2,out);
        well.ExecuteSimulation(nsubsteps,out);
        well.PostProcess(0);
        
        cout << "Average vertical stress " << well.GetCurrentConfig()->AverageVerticalStress() << std::endl;
        
        well.SetFluidModel(EPenetrating);
        well.ExecuteSimulation(nsubsteps,out);
        well.PostProcess(0);
        cout << "Penetrating fluid Average vertical stress " << well.GetCurrentConfig()->AverageVerticalStress() << std::endl;
        
        
//         well.EvolveBothPressures(4, WellPressure*1.2,reservoirPressure);
//         cout << "Higher well pressure Average vertical stress " << well.GetCurrentConfig()->AverageVerticalStress() << std::endl;
//         
//         well.EvolveBothPressures(4,WellPressure,reservoirPressure*0.8);
//         cout << "Lower reservoir pressure Average vertical stress " << well.GetCurrentConfig()->AverageVerticalStress() << std::endl;
//         
        TPZBFileStream save;
        save.OpenWrite("wellbore1.bin");
        well.Write(save);
    }
    
//QUEBRA AKI
    {
        TPZBFileStream read;
        read.OpenRead("wellbore1.bin");
        well.Read(read);
    }
    
}


//Reprodução do bug do erick para Philippe
void Config6()
{
    InitializePZLOG();
    gRefDBase.InitializeUniformRefPattern(EOned);
    gRefDBase.InitializeUniformRefPattern(ETriangle);
    gRefDBase.InitializeUniformRefPattern(EQuadrilateral);
    std::cout << std::setprecision(15);
    
    TPZWellBoreAnalysis well;
    
    STATE biotcoef = 0.6666666;
    well.SetBiotCoefficient(biotcoef);
    
    REAL reservoirPressure=57.2;
    well.SetReservoirPressure(reservoirPressure);
    
    
    REAL SH,Sh,SV;
    Sh=-83.5;
    SH=-99.8;
    SV=-85.9;
    TPZManVector<STATE,3> confinementTotal(3,0.);
    confinementTotal[0] = Sh;
    confinementTotal[1] = SH;
    confinementTotal[2] = SV;
    REAL WellPressure = 70.0; //66.6 61.1 57.2
    well.SetConfinementTotalStresses(confinementTotal, WellPressure);
    
    
    REAL innerradius = 4.25*0.0254;
    REAL outerradius = 100.;
    well.SetInnerOuterRadius(innerradius, outerradius);
    
    
    std::string output = "Config6.vtk";
    const int nsubsteps = 1;
    well.SetVtkOutPutName(output);
    REAL poisson = 0.2;
    REAL elast = 30000.;
    REAL cohesion = 1.e8; // Very very big
    REAL Phi = 1.5533430342749532; // 89 degrees

    well.SetMohrCoulombParameters(poisson, elast, cohesion, Phi, Phi);
    well.GetCurrentConfig()->fModel = EElastic;
    
    int Startfrom=0;

    if(Startfrom == 0){
        int porder = 1;
        int nrad=20;
        int ncircle = 40;
        REAL delx = 0.5*innerradius*M_PI_2/ncircle;
        TPZManVector<int,2> numdiv(2);
        numdiv[0] = nrad;
        numdiv[1] = ncircle;
        well.SetMeshTopology(delx, numdiv);
        well.GetCurrentConfig()->fWellConfig = EVerticalWell;
        well.GetCurrentConfig()->CreateGeometricMesh();
        
        well.GetCurrentConfig()->CreateComputationalMesh(porder);
        std::cout << "Average vertical stress " << well.GetCurrentConfig()->AverageVerticalStress() << std::endl;
        well.GetCurrentConfig()->CreatePostProcessingMesh();
        //REAL farfieldwork = well.GetCurrentConfig()->ComputeFarFieldWork();
        well.PostProcess(0);
        
    }

    if(Startfrom ==0){
        
        int nsteps = 1;
        int numnewton = 80;
        well.GetCurrentConfig()->ModifyWellElementsToQuadratic();
        TPZWellBoreAnalysis::TConfig::gNumThreads = 8;
        well.ExecuteInitialSimulation(nsteps, numnewton);
        well.PostProcess(0);
        TPZBFileStream save;
        save.OpenWrite("wellbore0.bin");
        well.Write(save);
        
    }
    
}

void Config7()
{
    //EVertical
    //ENonPenetrating
    InitializePZLOG();
    gRefDBase.InitializeUniformRefPattern(EOned);
    gRefDBase.InitializeUniformRefPattern(ETriangle);
    gRefDBase.InitializeUniformRefPattern(EQuadrilateral);
    std::cout << std::setprecision(15);
    
    TPBrAcidFunc acidfunc;
    acidfunc.StandardParameters();
    acidfunc.CalculaDerivados();
    REAL innerradius = 4.25*0.0254;
    REAL outerradius = 3.;
    TPZVec<REAL> x(2,0.);
    for (REAL r=innerradius; r <= outerradius; r += (outerradius-innerradius)/15.) {
        x[0] = r;
        TPZManVector<REAL,3> func(1);
        acidfunc.Execute(x, func);
        std::cout << "r = " << r << " Elast " << func[0] << std::endl;
    }

    TPZWellBoreAnalysis well;
    
#ifndef USING_TBB
    if (NumberOfThreads.get_value()>=0) {
        TPZWellBoreAnalysis::TConfig::gNumThreads=NumberOfThreads.get_value();
    }
#else
    int number_tbb=NumberOfThreads.get_value();
    if(number_tbb<=0)number_tbb=1;
    task_scheduler_init init(number_tbb);
#endif
    
    
    STATE biotcoef = 0.659;
    well.SetBiotCoefficient(biotcoef);
    
    REAL reservoirPressure=57.2;
    well.SetReservoirPressure(reservoirPressure);
    
    
    REAL SH,Sh,SV;
    Sh=-83.5;
    SH=-99.8;
    SV=-85.9;
    TPZManVector<STATE,3> confinementTotal(3,0.);
    confinementTotal[0] = Sh;
    confinementTotal[1] = SH;
    confinementTotal[2] = SV;
    REAL WellPressure = 57.2; //66.6 61.1 57.2
    well.SetConfinementTotalStresses(confinementTotal, WellPressure);
    
    
    well.SetInnerOuterRadius(innerradius, outerradius);
    
    
    std::string output = "Config7.vtk";
    well.SetVtkOutPutName(output);
    
    
    REAL sqj2_refine=0.0001;
    int Startfrom=1;
    const int nsubsteps = 5;
    if (Startfrom == 0)
    {
        well.SetInnerOuterRadius(innerradius, outerradius);
        
        REAL poisson = 0.203;
        REAL elast = 29269.;
        REAL A = 152.54;
        REAL B = 0.0015489;
        REAL C = 146.29;
        REAL R = 0.91969;
        REAL D = 0.018768;
        REAL W = 0.006605;
        
        bool modelMC =false;
        
        if (modelMC)
        {
            REAL cohesion = 13.;
            REAL Phi = 0.52;
            well.SetMohrCoulombParameters(poisson, elast, cohesion, Phi, Phi);
            
        }
        else
        {
            well.SetSanderDiMaggioParameters(poisson, elast, A, B, C, R, D, W);
            
            
        }
        
        int porder = 2;
        int nrad=20;
        int ncircle = 40;
        REAL delx = 0.5*innerradius*M_PI_2/ncircle;
        TPZManVector<int,2> numdiv(2);
        numdiv[0] = nrad;
        numdiv[1] = ncircle;
        well.SetMeshTopology(delx, numdiv);
        
        
        well.GetCurrentConfig()->fWellConfig = EVerticalWell;
        
        
        well.GetCurrentConfig()->CreateGeometricMesh();
        
        
        
        well.GetCurrentConfig()->CreateComputationalMesh(porder);
        
        
        well.GetCurrentConfig()->CreatePostProcessingMesh();
        
        well.PostProcess(0);
        
    }
    if (Startfrom ==0)
    {
        
        std::stringstream out;
        out << "Applying 0th pass\n";
        int nsteps = 5;
        int numnewton = 90;
        well.GetCurrentConfig()->ModifyWellElementsToQuadratic();
        well.ExecuteInitialSimulation(nsteps, numnewton);
        well.PRefineElementAbove(sqj2_refine, 3,out);
        well.DivideElementsAbove(sqj2_refine,out);
        well.ExecuteSimulation(nsubsteps,out);
        well.PostProcess(0);
        TPZBFileStream save;
        save.OpenWrite("Config7-0.bin");
        well.Write(save);
        
    }
    
    if (Startfrom ==1)
    {
        TPZBFileStream read;
        read.OpenRead("Config7-0.bin");
        well.Read(read);
    }
    
    
    if (Startfrom <=1)
    {
        std::cout << "\n ------- 1 -------- "<<std::endl;
        std::stringstream out;
        out << "Applying 1st pass\n";
        
        
        well.GetCurrentConfig()->fAcidParameters.StandardParameters();
        well.GetCurrentConfig()->ActivateAcidification();
        
        well.ExecuteSimulation(1,out);
        
        well.PostProcess(1);
        
        TPZBFileStream save;
        save.OpenWrite("Config7-1.bin");
        well.Write(save);
        
        
    }
}

void Config8()
{
    //EVertical
    //ENonPenetrating
    InitializePZLOG();
    gRefDBase.InitializeUniformRefPattern(EOned);
    gRefDBase.InitializeUniformRefPattern(ETriangle);
    gRefDBase.InitializeUniformRefPattern(EQuadrilateral);
    std::cout << std::setprecision(15);
    
    TPZWellBoreAnalysis well;
    
#ifndef USING_TBB
    if (NumberOfThreads.get_value()>=0) {
        TPZWellBoreAnalysis::TConfig::gNumThreads=NumberOfThreads.get_value();
    }
#else
    int number_tbb=NumberOfThreads.get_value();
    if(number_tbb<=0)number_tbb=1;
    task_scheduler_init init(number_tbb);
#endif
    
    
    STATE biotcoef = 0.659;
    well.SetBiotCoefficient(biotcoef);
    
    REAL reservoirPressure=57.2;
    well.SetReservoirPressure(reservoirPressure);
    
    
    REAL SH,Sh,SV;
    Sh=-83.5;
    SH=-99.8;
    SV=-85.9;
//    Sh=-57.2;
//    SH=-57.2;
//    SV=-57.2;
    TPZManVector<STATE,3> confinementTotal(3,0.);
    confinementTotal[0] = Sh;
    confinementTotal[1] = SH;
    confinementTotal[2] = SV;
    REAL WellPressure = 57.2; //66.6 61.1 57.2
    well.SetConfinementTotalStresses(confinementTotal, WellPressure);
    
    
    REAL innerradius = 4.25*0.0254;
    REAL outerradius = 3.;
    well.SetInnerOuterRadius(innerradius, outerradius);
    
    
    std::string output = "Config8.vtk";
    well.SetVtkOutPutName(output);
    
    
    REAL sqj2_refine=0.0001;
    int Startfrom=0;
    const int nsubsteps = 5;
    if (Startfrom == 0)
    {
        well.SetInnerOuterRadius(innerradius, outerradius);
        
        REAL poisson = 0.203;
        REAL elast = 29269.;
        
        well.SetElasticParameters(poisson, elast);
        
        int porder = 2;
        int nrad=20;
        int ncircle = 40;
        REAL delx = 0.5*innerradius*M_PI_2/ncircle;
        TPZManVector<int,2> numdiv(2);
        numdiv[0] = nrad;
        numdiv[1] = ncircle;
        well.SetMeshTopology(delx, numdiv);
        
        
        well.GetCurrentConfig()->fWellConfig = EVerticalWell;
        
        
        well.GetCurrentConfig()->CreateGeometricMesh();
        
        
        well.GetCurrentConfig()->CreateComputationalMesh(porder);
        
        {
            std::stringstream sout;
            well.PrintInitialConfiguration(sout);
            well.AppendExecutionLog(sout);
        }

        
        well.GetCurrentConfig()->CreatePostProcessingMesh();
        
        well.PostProcess(0);
        
    }
    if (Startfrom ==0)
    {
        
        int nsteps = 5;
        int numnewton = 90;
        well.ExecuteInitialSimulation(nsteps, numnewton);
        well.PostProcess(0);
        TPZBFileStream save;
        save.OpenWrite("Config8-0.bin");
        well.Write(save);
        
    }
    
    if (Startfrom ==1)
    {
        TPZBFileStream read;
        read.OpenRead("Config8-0.bin");
        well.Read(read);
    }
    
    
    if (Startfrom <=1)
    {
        std::cout << "\n ------- 1 -------- "<<std::endl;
        std::stringstream out;
        out << "Applying 1st pass\n";
        
        well.PRefineElementAbove(sqj2_refine, 3,out);
        well.DivideElementsAbove(sqj2_refine,out);
        well.ExecuteSimulation(nsubsteps,out);
        well.PostProcess(0);
        TPZBFileStream save;
        save.OpenWrite("Config1-1.bin");
        well.Write(save);
        
        well.AppendExecutionLog(out);
    }
    
    well.PrintExecutionLog(std::cout);
    
}

void Config9()
{
    //EVertical
    //ENonPenetrating
    InitializePZLOG();
    gRefDBase.InitializeUniformRefPattern(EOned);
    gRefDBase.InitializeUniformRefPattern(ETriangle);
    gRefDBase.InitializeUniformRefPattern(EQuadrilateral);
    std::cout << std::setprecision(15);
    
    TPZWellBoreAnalysis well;
    
#ifndef USING_TBB
    if (NumberOfThreads.get_value()>=0) {
        TPZWellBoreAnalysis::TConfig::gNumThreads=NumberOfThreads.get_value();
    }
#else
    int number_tbb=NumberOfThreads.get_value();
    if(number_tbb<=0)number_tbb=1;
    task_scheduler_init init(number_tbb);
#endif
    
    
    STATE biotcoef = 0.659;
    well.SetBiotCoefficient(biotcoef);
    
    REAL reservoirPressure=57.2;
    well.SetReservoirPressure(reservoirPressure);
    
    
    REAL SH,Sh,SV;
    Sh=-83.5;
    SH=-99.8;
    SV=-85.9;
    TPZManVector<STATE,3> confinementTotal(3,0.);
    confinementTotal[0] = Sh;
    confinementTotal[1] = SH;
    confinementTotal[2] = SV;
    REAL WellPressure = 57.2; //66.6 61.1 57.2
    well.SetConfinementTotalStresses(confinementTotal, WellPressure);
    
    
    REAL innerradius = 4.25*0.0254;
    REAL outerradius = 3.;
    well.SetInnerOuterRadius(innerradius, outerradius);
    

    std::string output = "Config9.vtk";
    well.SetVtkOutPutName(output);
    
    
    int Startfrom=0;
    if (Startfrom == 0)
    {
        well.SetInnerOuterRadius(innerradius, outerradius);
        
        REAL poisson = 0.203;
        REAL elast = 29269.;
        REAL A = 152.54;
        REAL B = 0.0015489;
        REAL C = 146.29;
        REAL R = 0.91969;
        REAL D = 0.018768;
        REAL W = 0.006605;
        
        REAL ESteel = 210000.;
        REAL poissonSteel = 0.30;
        
        
        bool modelMC =false;
        
        if (modelMC)
        {
            REAL cohesion = 13.;
            REAL Phi = 0.52;
            well.SetMohrCoulombParameters(poisson, elast, cohesion, Phi, Phi);
            
        }
        else
        {
            well.SetSanderDiMaggioParameters(poisson, elast, A, B, C, R, D, W);
            
            
        }
        
        int porder = 2;
        int nrad=20;
        int ncircle = 40;
        REAL delx = 0.5*innerradius*M_PI_2/ncircle;
        TPZManVector<int,2> numdiv(2);
        numdiv[0] = nrad;
        numdiv[1] = ncircle;
        well.SetMeshTopology(delx, numdiv);
        
        
        well.GetCurrentConfig()->fWellConfig = EVerticalWell;
        
        
        well.GetCurrentConfig()->CreateGeometricMesh();
        
        
        well.GetCurrentConfig()->CreateComputationalMesh(porder);
        
        
        well.GetCurrentConfig()->CreatePostProcessingMesh();
        
        well.PostProcess(0);
        
    }
    if (Startfrom ==0)
    {
        
        int nsteps = 5;
        int numnewton = 90;
        well.ExecuteInitialSimulation(nsteps, numnewton);
        
        well.PostProcess(0);
        
        TPZBFileStream save;
        save.OpenWrite("Config9-0.bin");
        well.Write(save);
        
    }
    
    well.PrintExecutionLog(std::cout);
}

#include "pzgeoelbc.h"
#include "tpzarc3d.h"

#include "TPZRefPattern.h"
#include "tpzgeoelrefpattern.h"
#include "tpzgeoblend.h"
TPZGeoMesh *CreatGeoMesh()
{

	TPZGeoMesh *geomesh = new TPZGeoMesh();

	geomesh->NodeVec().Resize(6);
	TPZVec<REAL> coord(2);
	coord[0] = 100.;
	coord[1] = 0.;
	geomesh->NodeVec()[0] = TPZGeoNode(0, coord, *geomesh);
	coord[0] = 200.;
	coord[1] = 0.;
	geomesh->NodeVec()[1] = TPZGeoNode(1, coord, *geomesh);
	coord[0] = 141.421;
	coord[1] = 141.421;
	geomesh->NodeVec()[2] = TPZGeoNode(2, coord, *geomesh);
	coord[0] = 0.;
	coord[1] = 200.;
	geomesh->NodeVec()[3] = TPZGeoNode(3, coord, *geomesh);

	coord[0] = 0.;
	coord[1] = 100.;
	geomesh->NodeVec()[4] = TPZGeoNode(4, coord, *geomesh);

	coord[0] = 70.7107;
	coord[1] = 70.7107;
	geomesh->NodeVec()[5] = TPZGeoNode(5, coord, *geomesh);

	//TPZGeoEl *gel[2];

	//TPZVec <long> TopoQuad(4);
	//TopoQuad[0] = 0;
	//TopoQuad[1] = 1;
	//TopoQuad[2] = 2;
	//TopoQuad[3] = 5;
	//long index;
	//gel[0] = geomesh->CreateGeoBlendElement(EQuadrilateral, TopoQuad, 1, index);

	//TopoQuad[0] = 2;
	//TopoQuad[1] = 3;
	//TopoQuad[2] = 4;
	//TopoQuad[3] = 5;
	//gel[1] = geomesh->CreateGeoBlendElement(EQuadrilateral, TopoQuad, 1, index);

	TPZGeoEl *gel[3];

	TPZVec <long> TopoQuad(4);
	TopoQuad[0] = 0;
	TopoQuad[1] = 1;
	TopoQuad[2] = 3;
	TopoQuad[3] = 4;
	long index;

	gel[0] = new TPZGeoElRefPattern< pzgeom::TPZGeoBlend< pzgeom::TPZGeoQuad> >(index, TopoQuad, 1, *geomesh);

	TPZVec <long> TopoArc(3);
	TopoArc[0] = 1;
	TopoArc[1] = 3;
	TopoArc[2] = 2;
	gel[1] =  new TPZGeoElRefPattern< pzgeom::TPZArc3D >(index, TopoArc, 1, *geomesh);

	TopoArc[0] = 4;
	TopoArc[1] = 1;
	TopoArc[2] = 5;
	gel[2] = new TPZGeoElRefPattern< pzgeom::TPZArc3D >(index, TopoArc, 1, *geomesh);


	geomesh->BuildConnectivity();

//	TPZGeoElBC t1(gel[0], 4, -1);//bottom

///	TPZGeoElBC t2(gel[0], 7, -2);//inner

//	TPZGeoElBC t3(gel[1], 5, -3);//top

	//TPZGeoElBC t4(gel[1], 6, -2);//inner

	geomesh->Print(std::cout);

	return geomesh;
}




#include "TPZRefPatternDataBase.h"
#include "TPZRefPatternTools.h"

#include <iostream>
#include <string>
#include <math.h>

#include "pzelasAXImat.h" 
#include "pzfstrmatrix.h"
#include "pzbstrmatrix.h"
#include "TPZParFrontStructMatrix.h"
#include "TPZInterfaceEl.h"


#include "pzlog.h"

#include "pzvec.h"
#include "pzstack.h"
#include "pzfmatrix.h"
#include "pzfstrmatrix.h"
#include "pzlog.h"

#include "pzgmesh.h"
#include "pzcmesh.h"
#include "pzcompel.h"
#include "pzgeoelside.h"
#include "TPZGeoLinear.h"
#include "pzgeopoint.h"
#include "tpzgeoblend.h"

#include <pzgengrid.h>
#include "TPZMatElasticity2D.h"
#include "TPZInterfaceEl.h"
#include "pzdiscgal.h"

#include "TPZRefPattern.h"
#include "tpzgeoelrefpattern.h"
#include "tpzcompmeshreferred.h"
#include "tpzautopointer.h"
#include "pzbndcond.h"
#include "pzanalysis.h"
#include <tpzarc3d.h>

#include "TPZParSkylineStructMatrix.h"
#include "pzstepsolver.h"
#include "pzstrmatrix.h"
#include "TPZFrontNonSym.h"
#include "TPZFrontSym.h"
#include "TPBSpStructMatrix.h"
#include "TPZSpStructMatrix.h"
#include "pzbstrmatrix.h"
#include "pzl2projection.h"

#include "pzpoisson3d.h"
#include "pzpoisson3dreferred.h"

#include "pzelasmat.h"
#include "pzmultiphysicselement.h"
#include "pzmultiphysicscompel.h"
#include "pzbuildmultiphysicsmesh.h"
#include "TPZSpStructMatrix.h"
#include "pzlog.h"
#include <iostream>
#include <string>
#include "TPZVTKGeoMesh.h"
#include "pzfunction.h"
#include "TPZReadGIDGrid.h"
#include "pzmultiphysicselement.h"
#include "TPZMultiphysicsInterfaceEl.h"
#include "TPZRefPatternDataBase.h"
#include "TPZRefPatternTools.h"

#include <cmath>
#include <set>
TPZGeoMesh * ShellGeoGen()
{
	int Qnodes = 6;

	TPZGeoMesh * gMesh = new TPZGeoMesh;

	//gMesh->InitializeRefPatterns();
	//gMesh->ImportRefPattern();

	gRefDBase.InitializeAllUniformRefPatterns();
	gRefDBase.InitializeRefPatterns();

	gMesh->NodeVec().Resize(Qnodes);

	TPZVec <long> TopolQuad(4);
	TPZVec <long> TopolArc(3);
	TPZVec <long> TopolLine(2);

	TPZVec<TPZGeoNode> Node(Qnodes);

	long id = 0;

	STATE a = 100.;
	STATE b = 200.;

	Node[id].SetNodeId(id);
	Node[id].SetCoord(0, a);//coord X
	Node[id].SetCoord(1, 0);//coord Y
	gMesh->NodeVec()[id] = Node[id];
	id++;

	Node[id].SetNodeId(id);
	Node[id].SetCoord(0, b);//coord X
	Node[id].SetCoord(1, 0);//coord Y
	gMesh->NodeVec()[id] = Node[id];
	id++;

	Node[id].SetNodeId(id);
	Node[id].SetCoord(0, 0.);//coord X
	Node[id].SetCoord(1, b);//coord Y
	gMesh->NodeVec()[id] = Node[id];
	id++;

	Node[id].SetNodeId(id);
	Node[id].SetCoord(0, 0.);//coord X
	Node[id].SetCoord(1, a);//coord Y
	gMesh->NodeVec()[id] = Node[id];
	id++;

	Node[id].SetNodeId(id);
	Node[id].SetCoord(0, a*sqrt(2.) / 2.);//coord X
	Node[id].SetCoord(1, a*sqrt(2.) / 2.);//coord Y
	gMesh->NodeVec()[id] = Node[id];
	id++;

	Node[id].SetNodeId(id);
	Node[id].SetCoord(0, b*sqrt(2.) / 2.);//coord X
	Node[id].SetCoord(1, b*sqrt(2.) / 2.);//coord Y
	gMesh->NodeVec()[id] = Node[id];
	id++;


	id = 0;
	TopolQuad[0] = 0;	TopolQuad[1] = 1;
	TopolQuad[2] = 2; 	TopolQuad[3] = 3;
	new TPZGeoElRefPattern< pzgeom::TPZGeoBlend< pzgeom::TPZGeoQuad> >(id, TopolQuad, 1, *gMesh);
	id++;

	TopolArc[0] = 0;
	TopolArc[1] = 3;
	TopolArc[2] = 4;
	new TPZGeoElRefPattern< pzgeom::TPZArc3D >(id, TopolArc, 1, *gMesh);
	id++;

	TopolArc[0] = 1;
	TopolArc[1] = 2;
	TopolArc[2] = 5;
	new TPZGeoElRefPattern< pzgeom::TPZArc3D >(id, TopolArc, 1, *gMesh);
	id++;

	TopolLine[0] = 0;
	TopolLine[1] = 1;
	new TPZGeoElRefPattern< pzgeom::TPZGeoLinear >(id, TopolLine, 1, *gMesh);
	id++;

	TopolLine[0] = 3;
	TopolLine[1] = 2;
	new TPZGeoElRefPattern< pzgeom::TPZGeoLinear >(id, TopolLine, 1, *gMesh);
	id++;

	gMesh->BuildConnectivity();
	set<int> SETmatDirectRef;
	int nh = 0;
	RefinamentoUniforme(*gMesh,nh);

	//for (int i = 0; i < 3; i++)
	//{
	//	int nel = gMesh->NElements();
	//	for (int iref = 0; iref < nel; iref++)
	//	{
	//		TPZVec<TPZGeoEl*> filhos;
	//		TPZGeoEl * gelP1 = gMesh->ElementVec()[iref];
	//		if (!gelP1) continue;
	//		//set<int> SETmatDirectRef;
	//		SETmatDirectRef.insert(1);
	//		//TPZRefPattern::RefineDirectional(gelP1, SETmatDirectRef);
	//		TPZRefPatternTools::RefineDirectional(gelP1, SETmatDirectRef);
	//	}
	//}
	
	ofstream arg("gmesh.txt");
	gMesh->Print(arg);

	return gMesh;
}
void RefinamentoUniforme(TPZGeoMesh & gMesh, int &nh)
{
	//int h =nh;
	for (int ref = 0; ref < nh; ref++)
	{// h indica o numero de refinamentos
		TPZVec<TPZGeoEl *> filhos;
		long n = gMesh.NElements();
		for (long i = 0; i < n; i++)
		{
			TPZGeoEl * gel = gMesh.ElementVec()[i];
			//int ind =  gel->MaterialId();
			gel->Divide(filhos);
		}//for i
	}//ref
}
/** FORCING FUNCTION */
void ForcingFunction2(const TPZVec<REAL> &x, TPZVec<STATE> &sol, TPZFMatrix<STATE> &dsol)
{
	std::cout << "x = " << x[0] << std::endl;
	std::cout << "y = " << x[1] << std::endl;
	//sol[1] = 100.*sin(M_PI*x[0] / 16.);
	//sol[1] = exp(-x[0] + 8.) / 10.;
	sol[1] = exp((-x[0] + 1.) / 0.1) / 10.;
	//sol[1] = 10000.;
	dsol.Zero();
}



TPZCompMesh *CreateCompMesh(TPZGeoMesh *geomesh)
{
	// Criacao da malha computacional
	TPZCompMesh *comp = new TPZCompMesh(geomesh);
	comp->SetDefaultOrder(2);
	TPZElastoPlasticAnalysis::SetAllCreateFunctionsWithMem(comp);

	int materialid = 1;
	bool planestrain = true;
	
	REAL elast = 210.;
	REAL poisson = 0.3;
	REAL K = elast / (3. * (1. -2.* poisson));
	REAL G = elast / (2. * (1. + poisson));
	REAL sigy = 0.24;
	TPZElasticResponse ER;
	ER.SetUp(elast, poisson);
	TPZPlasticStepPV<TPZYCVonMises, TPZElasticResponse> VM;
	VM.fYC.SetUp(K,G,sigy);
	VM.fER.SetUp(elast, poisson);
	//void TPZYCVonMises::ProjectSigmaDep(const TPZVec<STATE> &sigtrial, STATE kprev, TPZVec<STATE> &sigproj, STATE &kproj, TPZFMatrix<STATE> &GradSigma) const
	//TPZVec<REAL> pttrial(3, 0.),ptproj;
	//TPZFMatrix<STATE> Dep;
	//REAL sigyn1;
	//pttrial[0] = 13.1582;
	//pttrial[1] = 6.02444;
	//pttrial[2] = 5.75481;
	//VM.fYC.ProjectSigmaDep(pttrial,sigy,ptproj,sigyn1, Dep);
	//std::cout << ptproj[0] << std::endl;
	//std::cout << ptproj[1] << std::endl;
	//std::cout << ptproj[2] << std::endl;

	//TPZTensor<REAL> epst,sigma;
	//epst.XX() = 0.0025;
	//epst.YY() = 0.045;
	//epst.XY() = -0.012;
	//TPZFMatrix<REAL> dep(6, 6, 0.);
	//VM.ApplyStrainComputeDep(epst, sigma, dep);
	

	TPZMatElastoPlasticSest2D<TPZPlasticStepPV<TPZYCVonMises, TPZElasticResponse> > *PlasticVM = new TPZMatElastoPlasticSest2D< TPZPlasticStepPV<TPZYCVonMises, TPZElasticResponse> >(materialid, planestrain);
	PlasticVM->SetPlasticity(VM);
	comp->InsertMaterialObject(PlasticVM);

	// Creating four boundary condition
	TPZFMatrix<STATE> val1(2, 2, 0.), val2(2, 1, 0.);
	TPZMaterial *bcBottom, *bcRight, *bcTop;

	val2(1, 0) = 1.;
	bcBottom = PlasticVM->CreateBC(PlasticVM, -1, 3, val1, val2);

	val2(0, 0) = 1.;
	val2(1, 0) = 0.;
	bcRight = PlasticVM->CreateBC(PlasticVM, -3, 3, val1, val2);

	val2(1, 0) = 0.19209;
	val2(0, 0) = 0.;
	bcTop = PlasticVM->CreateBC(PlasticVM, -2, 5, val1, val2);

	//bcTop->SetForcingFunction(new TPZDummyFunction<STATE>(ForcingFunction2));

	comp->InsertMaterialObject(PlasticVM);
	comp->InsertMaterialObject(bcBottom);
	comp->InsertMaterialObject(bcRight);
	comp->InsertMaterialObject(bcTop);

	// Ajuste da estrutura de dados computacional
	comp->AutoBuild();
	//  comp->Print(cout);
	comp->AdjustBoundaryElements();
	//  comp->Print(cout);
	comp->CleanUpUnconnectedNodes();

	comp->SetName("Malha Computacional Original");

	return comp;
}

void ExecuteSimulation(TPZCompMesh cmesh, std::ostream &out);
#include "TPZVTKGeoMesh.h"
int main(int argc, char **argv)
{

	//TPZGeoMesh *gmesh = CreatGeoMesh();
	TPZGeoMesh *gmesh = ShellGeoGen();
	std::string vtkFile = "saida-tubo.vtk";
	std::ofstream argument("GeometicMeshsdf.txt");
	gmesh->Print(argument);
	std::ofstream Dummyfile(vtkFile);
	TPZVTKGeoMesh::PrintGMeshVTK(gmesh, Dummyfile, true);

	TPZCompMesh * cmesh = CreateCompMesh(gmesh);
	std::ofstream out("saida.txt");




	TPZStack<std::string> scalNames, vecNames;
	TPZWellBoreAnalysis::PostProcessVariables(scalNames, vecNames);



	ExecuteSimulation(*cmesh, out);


	std::cout << "HELLO" << std::endl;
	system("PAUSE");
    return 0;
}



void ExecuteSimulation(TPZCompMesh cmesh, std::ostream &out)
{

	int BCId = -3;
	TPZMaterial * mat = cmesh.FindMaterial(BCId);
	TPZBndCond * pBC = dynamic_cast<TPZBndCond *>(mat);

	TPZElastoPlasticAnalysis analysis(&cmesh, std::cout);

	TPZSkylineStructMatrix full(&cmesh);

	analysis.IdentifyEquationsToZero();

	TPZStepSolver<REAL> step;
	step.SetDirect(ELDLt);

	long neq = cmesh.NEquations();
	TPZVec<long> activeEquations;
	analysis.GetActiveEquations(activeEquations);
	TPZEquationFilter filter(neq);
	filter.SetActiveEquations(activeEquations);
	full.EquationFilter() = filter;
	analysis.SetStructuralMatrix(full);

	//step.SetDirect(ECholesky);
	analysis.SetSolver(step);

	int NumIter = 50;
	bool linesearch = true;
	bool checkconv = false;
	REAL tol = 1.e-5;
	bool conv;

    analysis.IterativeProcess(out, tol, NumIter, linesearch, checkconv, conv);
	TPZFMatrix<STATE> &sol = analysis.Mesh()->Solution();
	analysis.AcceptSolution();

}