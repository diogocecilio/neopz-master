#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include <QObject>
#include <QTemporaryFile>

//UI headers
#include "evaluatedw.h"
#include "parametersdw.h"
#include "loadparamsdw.h"
#include "geometryparamsdw.h"
#include "identifyellipsdw.h"
#include "pwbdw.h"
#include "initialconfigdialog.h"
#include "historydw.h"
#include "prefinedw.h"
#include "hrefinedw.h"
#include "fluiddw.h"
#include "poromechdw.h"
#include "meshdw.h"
#include "acidjobdw.h"
#include "casingdw.h"
#include "solutiondw.h"

//PZ headers
#include "WellBoreAnalysis.h"

//VTK headers
#include <vtkSmartPointer.h>
#include <vtkUnstructuredGrid.h>
#include <vtkActor.h>
#include <vtkDataSetMapper.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkContourFilter.h>
#include <vtkScalarBarActor.h>

namespace Ui {
class MainWindow;
}

class MainWindow : public QMainWindow
{
    Q_OBJECT
    
public:
    explicit MainWindow(QWidget *parent = 0);
    ~MainWindow();

    void openBin();
    void loadDefaultData();
    void dumpLogToWellboreObj();
    void loadLogFromWellboreObj();
    void reloadHistory();
    int writeHistory2(QString logMsg, bool append=false, int append_row=0);
    void composeScene();
    void initialZoomVTK();

private slots:
    void pwb_RunBtn();
    void historyClicked(int row);
    void changedMinMax();
    void postProcessVarChanged ();
    void initialConfigChanged();
    void redrawEllips();
    void redrawJ2();
    void redrawPolyChain();
    void redrawInnerRad();
    void BreakApplyBtn_clicked();
    void divideBtn_clicked();
    void prefineBtn_clicked();
    void calculateABBtn_clicked();
    void fluidModelChanged (int fluidModel);
    void ShowHistoryContextMenu(const QPoint &pos);
    void populate_postPVar2CBox(int idx);
    void changedTabs_updateCheckboxes (int idx);
    void changedTabs_updateCheckboxes (bool sts);
    void populate_GUIFields();
    void save_GUIFields();

    void on_gridCheck_toggled(bool checked);
    void on_solutionCheck_toggled(bool checked);
    void ellipsCheck_toggled(bool checked);
    void sqrtJ2Check_toggled(bool checked);
    void on_innerRadCheck_toggled(bool checked);
    void contourCheck_toggled(bool checked);
    void on_showFullCheck_toggled(bool checked);
    void polyChainCheck_toggled(bool checked);
    void on_scaleCheck_toggled(bool checked);
    void on_resetZoomBtn_clicked();
    void on_actionInitialConfig_triggered(bool checked);
    void on_actionRestart_triggered(bool checked);

    void on_actionOpen_triggered();

    void on_actionSave_triggered();

private:
    Ui::MainWindow *ui;

    // generate a discontinuous representation of the geometric grid
    vtkSmartPointer<vtkUnstructuredGrid> GenerateVTKGridDisc(double setZ = 0);

    // generate the data of any variable
    void computeVisualPostProcess(TPZVec<STATE> &vecvar, QString varName, QString varType);

    //
    void addPostProcVar(int i);

    //Docks
    evaluateDW *evaluateDock;
    parametersDW *parametersDock;
    loadparamsDW *loadparamsDock;
    geometryparamsDW *geometryparamsDock;
    identifyellipsDW *identifyellipsDock;
    pwbDW *pwbDock;
    historyDW *historyDock;
    prefineDW *prefineDock;
    hrefineDW *hrefineDock;
    poromechDW *poromechDock;
    fluidDW *fluidDock;
    meshDW *meshDock;
    acidjobDW *acidjobDock;
    casingDW *casingDock;
    solutionDW *solutionDock;

    //Dialogs
    initialConfigDialog *initialDialog;

    //VTK object pointers
    vtkSmartPointer<vtkUnstructuredGrid> aGrid;

    vtkSmartPointer<vtkDataSetMapper> mapperMesh;
    vtkSmartPointer<vtkDataSetMapper> mapperWire;
    vtkSmartPointer<vtkPolyDataMapper> mapperEllips;
    vtkSmartPointer<vtkPolyDataMapper> mapperPolychain;
    vtkSmartPointer<vtkPolyDataMapper> contMapper;
    vtkSmartPointer<vtkPolyDataMapper> contMapper2;
    vtkSmartPointer<vtkPolyDataMapper> mapperInnerRad;

    vtkSmartPointer<vtkActor> actorMesh;
    vtkSmartPointer<vtkActor> actorWire;
    vtkSmartPointer<vtkActor> actorEllips;
    vtkSmartPointer<vtkActor> actorPolychain;
    vtkSmartPointer<vtkActor> actorCont;
    vtkSmartPointer<vtkActor> actorCont2;
    vtkSmartPointer<vtkActor> actorInnerRad;

    vtkSmartPointer<vtkRenderer> renderer;
    vtkSmartPointer<vtkRenderWindow> renderer_window;

    vtkSmartPointer<vtkRenderWindowInteractor> interactor;

    vtkSmartPointer<vtkContourFilter> contours;
    vtkSmartPointer<vtkContourFilter> contours2;

    vtkSmartPointer<vtkScalarBarActor> scalarBar;

    //WellboreAnalysis
    TPZWellBoreAnalysis fwellb;
    TPZWellBoreAnalysis::TConfig *fLoadedConfig;

    //Units
    QTemporaryFile fTmpDefsUnitsFile;

    //Min / Max values used into scaleBar
    double minVal;
    double maxVal;
};

#endif // MAINWINDOW_H
