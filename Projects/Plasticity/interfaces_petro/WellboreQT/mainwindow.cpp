#include "mainwindow.h"
#include "ui_mainwindow.h"

#include "ui_parametersdw.h"
#include "ui_evaluatedw.h"
#include "ui_geometryparamsdw.h"
#include "ui_identifyellipsdw.h"
#include "ui_loadparamsdw.h"
#include "ui_pwbdw.h"
#include "ui_initialconfigdialog.h"
#include "ui_historydw.h"
#include "ui_prefinedw.h"
#include "ui_hrefinedw.h"
#include "ui_poromechdw.h"
#include "ui_fluiddw.h"
#include "ui_meshdw.h"
#include "ui_acidjobdw.h"
#include "ui_casingdw.h"
#include "ui_solutiondw.h"

#include "wellboreinteractorstyle.h"

//VTK headers
#include <vtkVersion.h>
#include <vtkScalarBarWidget.h>
#include <vtkLookupTable.h>
#include <vtkScalarBarActor.h>
#include <vtkSmartPointer.h>
#include <vtkProperty.h>
#include <vtkContourWidget.h>
#include <vtkOrientedGlyphContourRepresentation.h>
#include <vtkRenderer.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkCommand.h>
//#include <vtkRegressionTestImage.h>
#include <vtkDebugLeaks.h>
//#include <vtkTestUtilities.h>
#include <vtkCamera.h>
#include <vtkPlane.h>
#include <vtkPolyData.h>
#include <vtkCellArray.h>
#include <vtkPoints.h>
#include <vtkMath.h>
#include <vtkWidgetEvent.h>
#include <vtkWidgetEventTranslator.h>
#include <vtkContourFilter.h>
#include <vtkImageData.h>
#include <vtkGenericDataObjectReader.h>
#include <vtkStructuredGrid.h>
#include <vtkSmartPointer.h>
#include <vtkPolyData.h>
#include <string>
#include <vtkAlgorithmOutput.h>
#include <vtkDataSetReader.h>
#include <vtkDataSet.h>
#include <vtkDataObject.h>
#include <vtkPointData.h>
#include <vtkCellData.h>
#include <vtkDataArray.h>
#include <vtkIndent.h>
#include <vtkPolyDataMapper.h>
#include <vtkStructuredGridGeometryFilter.h>
#include <vtkUnstructuredGridGeometryFilter.h>
#include <vtkUnstructuredGrid.h>
#include <vtkUnstructuredGridReader.h>
#include <vtkDataSetMapper.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkRenderer.h>
#include <vtkProperty.h>
#include <vtkSmartPointer.h>
#include <vtkRendererCollection.h>
#include <vtkPointPicker.h>
#include <vtkSphereSource.h>
#include <vtkPolyDataMapper.h>
#include <vtkActor.h>
#include <vtkObjectFactory.h>
#include <vtkRenderWindow.h>
#include <vtkRenderer.h>
#include <vtkPolyData.h>
#include <vtkPointSource.h>
#include <vtkInteractorStyle.h>
#include <vtkInteractorStyleTrackballCamera.h>
#include <vtkInteractorStyleRubberBandZoom.h>
#include <vtkInteractorStyleTrackballCamera.h>
#include <vtkInteractorStyleSwitch.h>
#include <vtkInteractorStyleRubberBand2D.h>
#include <vtkInteractorStyleRubberBand3D.h>
#include <vtkInteractorStyleJoystickActor.h>
#include <vtkAreaPicker.h>
#include <vtkExtractGeometry.h>
#include <vtkDataSetMapper.h>
#include <vtkUnstructuredGrid.h>
#include <vtkVertexGlyphFilter.h>
#include <vtkIdFilter.h>
#include <vtkCoordinate.h>
#include <vtkViewport.h>
#include <vtkOrientationMarkerWidget.h>
#include <vtkAxesActor.h>
#include <vtkTransform.h>
#include <vtkQuad.h>
#include <vtkTriangle.h>
#include <vtkDataSetMapper.h>
#include <vtkCellArray.h>
#include <vtkDoubleArray.h>
#include <vtkCaptionActor2D.h>
#include <vtkTextProperty.h>
#include <vtkCamera.h>
#include <vtkActor2DCollection.h>

////////////////////////////////////////////////////////////////////////
#define DEBUG

//QT headers
//#ifdef DEBUG
#include <QDebug>
//#endif
#include <QFile>
#include <QMessageBox>
#include <QFileDialog>
#include <QTemporaryFile>
#include <QTextStream>

int startfrom=0;

//PZ headers
#include "TPZRefPatternDataBase.h"
#include "TPZSandlerDimaggio.h"
#include "TPZProjectEllipse.h"
#include "pzmaterial.h"
#include "pzlog.h"
#include "pzvec.h"
#include "pzgmesh.h"
#include "pzstack.h"
#include "pzfmatrix.h"
#include "pzmanvector.h"
#include "pzbfilestream.h"
#include "TPBrBiotForce.h"

//#define TOBIN

#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("interfaces_petro.mainwindow"));
#endif

MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::MainWindow)
{
    // Creating a temporary file with resource file definition.units contents.
    // Temporary file will be removed when program ends!
    QFile unitsFile (":/units/definitions.units");
    if (unitsFile.open(QIODevice::ReadOnly | QIODevice::Text)) {
        if (fTmpDefsUnitsFile.open()) {
            QTextStream in(&unitsFile);
            QString line = in.readLine();
            QTextStream out(&fTmpDefsUnitsFile);

            while (!line.isNull()) {
                out << line << endl;
                line = in.readLine();
            }
            fTmpDefsUnitsFile.close();
        }
        unitsFile.close();
        // saving its name into TPBrUnitInput
        TPBrUnitInput::fTmpDefsUnitsFile = fTmpDefsUnitsFile.fileName();
    }

    //Set up window components
    ui->setupUi(this);

    setWindowTitle("SEST2D");
    // Dockwidgets setup
//    evaluateDock = new evaluateDW(this);
//    evaluateDock->setWindowTitle("Evaluate");
//    addDockWidget(Qt::LeftDockWidgetArea, evaluateDock);

    // default min / max values
    minVal = 0; //ui->minBar->text().toDouble();
    maxVal = 1; //ui->maxBar->text().toDouble();

    identifyellipsDock = new identifyellipsDW(this);
    identifyellipsDock->setWindowTitle(tr("Breakout"));
    identifyellipsDock->setTitleBarWidget(new QWidget(this));
    addDockWidget(Qt::LeftDockWidgetArea, identifyellipsDock);

    prefineDock = new prefineDW(this);
    prefineDock->setWindowTitle(tr("P Ref"));
    prefineDock->setTitleBarWidget(new QWidget(this));
    addDockWidget(Qt::LeftDockWidgetArea, prefineDock);

    hrefineDock = new hrefineDW(this);
    hrefineDock->setWindowTitle(tr("H Ref"));
    hrefineDock->setTitleBarWidget(new QWidget(this));
    addDockWidget(Qt::LeftDockWidgetArea, hrefineDock);

    acidjobDock = new acidjobDW(this);
    acidjobDock->setWindowTitle(tr("Acid job"));
    acidjobDock->setTitleBarWidget(new QWidget(this));
    addDockWidget(Qt::LeftDockWidgetArea, acidjobDock);

    pwbDock = new pwbDW(this);
    pwbDock->setWindowTitle(tr("Pwb"));
    pwbDock->setTitleBarWidget(new QWidget(this));
    addDockWidget(Qt::LeftDockWidgetArea, pwbDock);

    casingDock = new casingDW(this);
    casingDock->setWindowTitle(tr("Casing"));
    casingDock->setTitleBarWidget(new QWidget(this));
    addDockWidget(Qt::LeftDockWidgetArea, casingDock);

    solutionDock = new solutionDW(this);
    solutionDock->setWindowTitle(tr("Solution"));
    solutionDock->setTitleBarWidget(new QWidget(this));
    addDockWidget(Qt::LeftDockWidgetArea, solutionDock);

    tabifyDockWidget(pwbDock, identifyellipsDock);
    tabifyDockWidget(identifyellipsDock, prefineDock);
    tabifyDockWidget(prefineDock, hrefineDock);
    tabifyDockWidget(hrefineDock, acidjobDock);
    tabifyDockWidget(acidjobDock, casingDock);
    tabifyDockWidget(casingDock, solutionDock);
    pwbDock->raise(); //selected tab

    historyDock = new historyDW(this);
    historyDock->setWindowTitle(tr("History"));
    historyDock->setTitleBarWidget(new QWidget(this));
    addDockWidget(Qt::LeftDockWidgetArea, historyDock);

    setTabPosition(Qt::LeftDockWidgetArea, QTabWidget::North);

    // Dialog
    initialDialog = new initialConfigDialog(this);
    initialDialog->setTabPosition(Qt::LeftDockWidgetArea, QTabWidget::North);

    parametersDock = new parametersDW(initialDialog);
    parametersDock->setWindowTitle(tr("Material"));
    parametersDock->setTitleBarWidget(new QWidget(this));
    initialDialog->addDockWidget(Qt::LeftDockWidgetArea, parametersDock);

    geometryparamsDock = new geometryparamsDW(initialDialog);
    geometryparamsDock->setWindowTitle(tr("Geometry"));
    geometryparamsDock->setTitleBarWidget(new QWidget(this));
    initialDialog->addDockWidget(Qt::LeftDockWidgetArea, geometryparamsDock);

    meshDock = new meshDW(initialDialog);
    meshDock->setWindowTitle(tr("Mesh"));
    meshDock->setTitleBarWidget(new QWidget(this));
    initialDialog->addDockWidget(Qt::LeftDockWidgetArea, meshDock);

    loadparamsDock = new loadparamsDW(initialDialog);
    loadparamsDock->setWindowTitle(tr("Load stress"));
    loadparamsDock->setTitleBarWidget(new QWidget(this));
    //loadparamsDock->setFixedHeight(250); //gambis: to show all tabs
    //loadparamsDock->setFixedWidth(550); //gambis: to show all tabs
    loadparamsDock->ui->label_pr0->setText(QString("Pr<SUB>0</SUB>"));
    initialDialog->addDockWidget(Qt::LeftDockWidgetArea, loadparamsDock);

    //fluidDock = new fluidDW(initialDialog);
    //fluidDock->setWindowTitle(tr("Fluid"));
    //initialDialog->addDockWidget(Qt::LeftDockWidgetArea, fluidDock);

    //poromechDock = new poromechDW(initialDialog);
    //poromechDock->setWindowTitle(tr("PoroMech"));
    //initialDialog->addDockWidget(Qt::LeftDockWidgetArea, poromechDock);

    initialDialog->tabifyDockWidget(loadparamsDock, geometryparamsDock);
    initialDialog->tabifyDockWidget(geometryparamsDock, meshDock);
    initialDialog->tabifyDockWidget(meshDock, parametersDock);
    //initialDialog->tabifyDockWidget(parametersDock, fluidDock);
    //initialDialog->tabifyDockWidget(fluidDock, poromechDock);
    initialDialog->setStyleSheet("QTabBar::tab { height: 30px; width: 85px; }"); //size of tabs
    loadparamsDock->raise(); //selected tab

    // Adding solution vars

    ui->postPVar1CBox->addItem("Strain","Strain");
    ui->postPVar1CBox->addItem("Elastic Strain","ElasticStrain");
    ui->postPVar1CBox->addItem("Plastic Strain","PlasticStrain");
    ui->postPVar1CBox->addItem("Displacement","Displacement");
    ui->postPVar1CBox->addItem("Simulation","Simulation");
    ui->postPVar1CBox->addItem("Total Stress","TotalStress");
    ui->postPVar1CBox->addItem("Pore Pressure","PorePressure");
    ui->postPVar1CBox->addItem("Effective Stress","EffectiveStress");
    ui->postPVar1CBox->addItem("Yield Surface","YieldSurface");
    ui->postPVar1CBox->addItem("Material","Material");

    populate_postPVar2CBox(0);

    // Populate second postprocess combo box
    connect (ui->postPVar1CBox, SIGNAL(currentIndexChanged(int)), this, SLOT(populate_postPVar2CBox(int)));
    connect (ui->postPVar2CBox, SIGNAL(currentIndexChanged(int)), this, SLOT(postProcessVarChanged()));



    //LOGOS
    //ui->logoPetro->setText("<img src=':imgs/LogoPetro.png' height='60' />");
    //ui->logoLabmec->setText("<img src=':imgs/LogoLabMec.png' height='40' />");

    // Connects
    connect (prefineDock->ui->prefineBtn, SIGNAL(clicked()), this, SLOT(prefineBtn_clicked()));
    connect (hrefineDock->ui->hrefineBtn, SIGNAL(clicked()), this, SLOT(divideBtn_clicked()));
    connect (identifyellipsDock->ui->BreakApplyBtn, SIGNAL(clicked()), this, SLOT(BreakApplyBtn_clicked()));
    connect (identifyellipsDock->ui->calculateABBtn, SIGNAL(clicked()), this, SLOT(calculateABBtn_clicked()));
    connect (identifyellipsDock->ui->fa, SIGNAL(editingFinished()), this, SLOT(redrawEllips()));
    connect (identifyellipsDock->ui->fb, SIGNAL(editingFinished()), this, SLOT(redrawEllips()));
    connect (identifyellipsDock->ui->sqrtJ2, SIGNAL(editingFinished()), this, SLOT(redrawJ2()));
    connect (identifyellipsDock->ui->sqrtJ2, SIGNAL(editingFinished()), this, SLOT(redrawPolyChain()));
    connect (pwbDock->ui->runBtn, SIGNAL(clicked()), this, SLOT(pwb_RunBtn()));
    connect (historyDock->ui->historyList, SIGNAL(currentRowChanged(int)), this, SLOT(historyClicked(int)));
    // connecting combo box signal
    // connecting enable RunButton signal
    connect (initialDialog->ui->okBtn, SIGNAL(clicked()), this, SLOT(initialConfigChanged()) );
    // connecting min and max scale events
    connect(ui->minBar, SIGNAL(editingFinished()), this, SLOT(changedMinMax()));
    connect(ui->maxBar, SIGNAL(editingFinished()), this, SLOT(changedMinMax()));
    connect(ui->autoBar, SIGNAL(toggled(bool)), this, SLOT(postProcessVarChanged()));
    connect(ui->autoBar, SIGNAL(toggled(bool)), ui->minBar, SLOT(setDisabled(bool)));
    connect(ui->autoBar, SIGNAL(toggled(bool)), ui->maxBar, SLOT(setDisabled(bool)));
    // connecting fluidModel Radiobuttons
    connect(pwbDock, SIGNAL(fluidModelChanged(int)), this, SLOT(fluidModelChanged(int)));


    //PZ
    gRefDBase.InitializeAllUniformRefPatterns();

    fLoadedConfig = fwellb.GetCurrentConfig();

    loadDefaultData();
    populate_GUIFields();

    composeScene();

    // Connect master combo box units to all unit boxes inside each dockwidget
    // Change master unit, change all unit boxes
    connect (initialDialog->ui->pressureCBox, SIGNAL(currentIndexChanged(int)), loadparamsDock->ui->SH, SLOT(changeUnit(int)));
    connect (initialDialog->ui->pressureCBox, SIGNAL(currentIndexChanged(int)), loadparamsDock->ui->Sh, SLOT(changeUnit(int)));
    connect (initialDialog->ui->pressureCBox, SIGNAL(currentIndexChanged(int)), loadparamsDock->ui->SV, SLOT(changeUnit(int)));
    connect (initialDialog->ui->pressureCBox, SIGNAL(currentIndexChanged(int)), loadparamsDock->ui->porePressure, SLOT(changeUnit(int)));
    connect (initialDialog->ui->pressureCBox, SIGNAL(currentIndexChanged(int)), loadparamsDock->ui->wellborePressure, SLOT(changeUnit(int)));
    connect (initialDialog->ui->dimensionCBox, SIGNAL(currentIndexChanged(int)), loadparamsDock->ui->TVD, SLOT(changeUnit(int)));
    connect (initialDialog->ui->dimensionCBox, SIGNAL(currentIndexChanged(int)), geometryparamsDock->ui->r_int, SLOT(changeUnit(int)));
    connect (initialDialog->ui->dimensionCBox, SIGNAL(currentIndexChanged(int)), geometryparamsDock->ui->r_ext, SLOT(changeUnit(int)));

    // CONNECT ALL COMPONENTS THAT MAY HOLD PPG UNIT TYPE
    connect (loadparamsDock->ui->TVD, SIGNAL(valueChanged(double)), loadparamsDock->ui->Sh, SLOT(reload(double)));
    connect (loadparamsDock->ui->TVD, SIGNAL(valueChanged(double)), loadparamsDock->ui->SH, SLOT(reload(double)));
    connect (loadparamsDock->ui->TVD, SIGNAL(valueChanged(double)), loadparamsDock->ui->SV, SLOT(reload(double)));
    connect (loadparamsDock->ui->TVD, SIGNAL(valueChanged(double)), loadparamsDock->ui->porePressure, SLOT(reload(double)));
    connect (loadparamsDock->ui->TVD, SIGNAL(valueChanged(double)), loadparamsDock->ui->wellborePressure, SLOT(reload(double)));

    connect (loadparamsDock->ui->TVD, SIGNAL(valueChanged(double)), pwbDock->ui->pwb, SLOT(reload(double)));
    connect (loadparamsDock->ui->TVD, SIGNAL(valueChanged(double)), pwbDock->ui->pres, SLOT(reload(double)));
    // force all connects to act and update its TVD values
    loadparamsDock->ui->TVD->setText("2000");

    // Habilita menu de contexto
    historyDock->ui->historyList->setContextMenuPolicy(Qt::CustomContextMenu);
    connect (historyDock->ui->historyList, SIGNAL(customContextMenuRequested(const QPoint&)), this, SLOT(ShowHistoryContextMenu(const QPoint&)));


    // Get list of tabs (dockwidgets) stacked together
    QList<QTabBar *> tabList = this->findChildren<QTabBar *>();
    QTabBar *tabBar = tabList.at(1); // Checar e talvez mudar quando adicionar novas tablists
    // connect it to a slot, to get when tabs selection changes
    connect (tabBar, SIGNAL(currentChanged(int)), this, SLOT(changedTabs_updateCheckboxes(int)));

    connect (identifyellipsDock->ui->ellipsCheck, SIGNAL(toggled(bool)), this, SLOT(changedTabs_updateCheckboxes(bool)));
    connect (identifyellipsDock->ui->polyChainCheck, SIGNAL(toggled(bool)), this, SLOT(changedTabs_updateCheckboxes(bool)));
    connect (identifyellipsDock->ui->contourCheck, SIGNAL(toggled(bool)), this, SLOT(changedTabs_updateCheckboxes(bool)));
    connect (prefineDock->ui->polyChainCheck, SIGNAL(toggled(bool)), this, SLOT(changedTabs_updateCheckboxes(bool)));
    connect (prefineDock->ui->contourCheck, SIGNAL(toggled(bool)), this, SLOT(changedTabs_updateCheckboxes(bool)));
    connect (hrefineDock->ui->polyChainCheck, SIGNAL(toggled(bool)), this, SLOT(changedTabs_updateCheckboxes(bool)));
    connect (hrefineDock->ui->contourCheck, SIGNAL(toggled(bool)), this, SLOT(changedTabs_updateCheckboxes(bool)));


    changedTabs_updateCheckboxes (true);

    pwbDock->ui->runBtn->setEnabled(false);
}

MainWindow::~MainWindow()
{
    delete ui;
}

void MainWindow::initialConfigChanged() {

    // TODO: VALIDATE ALL INITIAL DATA HERE

    int configCount = fwellb.GetConfigListSize();

    // existe uma simulacao em aberto, reinicia a simulacao e seta com novos campos do initial config
    if (configCount != 0) {

        // TODO: DIALOGO INFORMANDO QUE A SIMULACAO CORRENTE SERA PERDIDA

        TPZWellBoreAnalysis a;
        fwellb = a;

        fLoadedConfig = fwellb.GetCurrentConfig();

        ui->logList->clear();
    }

    save_GUIFields();

    fLoadedConfig = fwellb.GetCurrentConfig();

    composeScene();

    //Enable Run/Simulate button
    initialDialog->hide();

    pwbDock->ui->pwb->setText(QVariant(fwellb.GetCurrentConfig()->fWellborePressure).toString());
    pwbDock->ui->pres->setText(QVariant(fwellb.GetCurrentConfig()->fReservoirPressure).toString());

    pwb_RunBtn();
    pwbDock->ui->runBtn->setEnabled(true);

    initialZoomVTK();

    ui->gridCheck->setEnabled(true);
    ui->scaleCheck->setEnabled(true);
    ui->innerRadCheck->setEnabled(true);
    ui->solutionCheck->setEnabled(true);
    ui->gridCheck->setChecked(true);
    ui->scaleCheck->setChecked(true);
    ui->innerRadCheck->setChecked(true);
    ui->solutionCheck->setChecked(true);
}

int MainWindow::writeHistory2 (QString logMsg, bool append, int append_row) {

    int rowCount = ui->logList->count();
    int row = 0;
    if (append == true && append_row < rowCount) {
        QString new_log (ui->logList->item(append_row)->text());
        new_log.append(logMsg);
        ui->logList->item(append_row)->setText(new_log);
#ifdef DEBUG
        qDebug() << new_log;
#endif
        row = append_row;
    }
    else {
        ui->logList->addItem(logMsg);
#ifdef DEBUG
        qDebug() << logMsg;
#endif
        row = rowCount;
    }
    ui->logList->scrollToBottom();
    ui->logList->repaint();
    QApplication::processEvents();

    return row;
}

void MainWindow::reloadHistory () {
#ifdef DEBUG
        qDebug() << "Reload history start";
#endif
    historyDock->ui->historyList->clear();

    //TODO: MELHORAR A LOGICA TOSCA ABAIXO

    //popula a lista
    // STORE CONFIG INDEX AT DATA-ROLE=5 FROM HISTORY ITEM
    for (int i = 0; i < fwellb.GetConfigListSize(); i++) {
        QString itemTitle (fwellb.GetConfig(i)->fHistoryLog.c_str());
        QListWidgetItem *newItem = new QListWidgetItem (itemTitle);
        newItem->setData(5,QVariant(i));
        historyDock->ui->historyList->addItem(newItem);
#ifdef DEBUG
        qDebug() << ">>>>>>>>>>>>>>>>>>>>>>>>>>> " << historyDock->ui->historyList->item(i)->text()
                 << "  id " << historyDock->ui->historyList->item(i)->data(5);
#endif
    }
    historyDock->ui->historyList->addItem(tr("Current Config"));


    int idx_parent_step = 1;
    //Atualiza itens pais com id dos seus ultimos filhos
    for (int i = 2; i < historyDock->ui->historyList->count(); i++) {
        QString itemTitle (historyDock->ui->historyList->item(i)->text());
        //se nao tiver identacao eh item pai
        if (itemTitle.toStdString().c_str()[0] != ' ') {
#ifdef DEBUG
            qDebug() << "######### PARENT " << idx_parent_step << " aponta para LAST " << i-1;
#endif

            QListWidgetItem *parentItem = historyDock->ui->historyList->item(idx_parent_step);
            parentItem->setData(5,QVariant(i-1));

            idx_parent_step = i;
        }
    }

    //Atualiza apontador do ultimo item da lista (current config) aponta para ultima config
    QListWidgetItem *parentItem = historyDock->ui->historyList->item(historyDock->ui->historyList->count()-1);
    parentItem->setData(5,QVariant(fwellb.GetConfigListSize()-1));

#ifdef DEBUG
    for (int i = 0; i < historyDock->ui->historyList->count(); i++) {
        qDebug() << ">>>>>>>>>>NOVA>>>>>>>>>>>>> " << historyDock->ui->historyList->item(i)->text()
                 << "  id " << historyDock->ui->historyList->item(i)->data(5);
    }
#endif

    historyDock->ui->historyList->scrollToBottom();
    historyDock->ui->historyList->repaint();
    QApplication::processEvents();

#ifdef DEBUG
        qDebug() << "Reload history done";
#endif
}

void MainWindow::loadDefaultData() {

        //SET DEFAULT VALUES INTO INTERFACE's FIELDS

        double inner = 4.25*0.0254;
        double outer = 3.;
        fwellb.SetInnerOuterRadius(inner, outer);

        REAL biot = 0.659; //0.; NAO PODE SER MAIOR OU IGUAL A HUM!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        fwellb.SetBiotCoefficient(biot);

        TPZManVector<STATE,3> confinement(3,0.);
        confinement[0] = -83.5; //-45.9;//-44.3;// MPa
        confinement[1] = -99.8; //-62.1;//-58.2;
        confinement[2] = -85.9; //-48.2;//-53.8;
        double WellBorePressure = 57.2; // 19.5;//29.3; //PWB
        double PorePressure = 57.2;

        //for(int i=0;i<3;i++)confinement[i]-=(biot*confinement[i]);
        fwellb.SetConfinementTotalStresses(confinement, WellBorePressure);
        fwellb.SetReservoirPressure(PorePressure);

        REAL poissonM = 0.203;
        REAL elastM = 29269.;
        REAL angleM = 0.52; //0.25;
        REAL cohesionM = 12.2; //7;
        fwellb.SetMohrCoulombParameters(poissonM, elastM, cohesionM, angleM, angleM);

        REAL poisson = 0.203;
        REAL elast = 29269.;
        REAL A = 152.54;
        REAL B = 0.0015489;
        REAL C = 146.29;
        REAL R = 0.91969;
        REAL D = 0.018768;
        REAL W = 0.006605;
        fwellb.SetSanderDiMaggioParameters(poisson, elast, A, B, C, R, D, W);

        int n_rad = 40;
        int n_circ = 20;
        meshDock->ui->n_radial->setText(QVariant(n_rad).toString());
        meshDock->ui->n_circ->setText(QVariant(n_circ).toString());

        REAL first_el_size = 0.5;
        meshDock->ui->first_el_size->setText(QVariant(first_el_size).toString());
        REAL Delx = first_el_size*inner*M_PI_2/n_circ;

        TPZManVector<int,2> numdiv(2);
        numdiv[0] = n_rad; //divisoes em X
        numdiv[1] = n_circ; //divisoes em Y
        fwellb.SetMeshTopology(Delx, numdiv);

        identifyellipsDock->ui->fa->setText(QVariant(inner*1.1).toString()); //("Not Computed");
        identifyellipsDock->ui->fb->setText(QVariant(inner*0.9).toString()); //("Not Computed");
        identifyellipsDock->ui->sqrtJ2->setText(QVariant(0.0007).toString());//("Not Computed");

        prefineDock->ui->p_order->setText("2");
        prefineDock->ui->isoValue->setText("0.0001");
        hrefineDock->ui->isoValue->setText("0.0001");

        loadparamsDock->ui->porePressure->setText(QVariant(PorePressure).toString());
        loadparamsDock->ui->wellborePressure->setText(QVariant(WellBorePressure).toString());

        //loadparamsDock->ui->num_steps->setText("3");

        //cleaning up history
        historyDock->ui->historyList->clear();
        ui->logList->clear();
//        pwbDock->ui->runBtn->setDisabled(true);

//        //Updating UI
//        ui->gridCheck->setChecked(true);
//        //ui->ellipsCheck->setChecked(true);
//        ui->solutionCheck->setChecked(true);
//        //ui->sqrtJ2Check->setChecked(true);
//        //ui->contourCheck->setChecked(true);
//        //ui->polyChainCheck->setChecked(true);
//        ui->scaleCheck->setChecked(true);
//        ui->innerRadCheck->setChecked(true);

//        //Disabling checks UI
//        ui->gridCheck->setEnabled(false);
//        //ui->ellipsCheck->setEnabled(false);
//        ui->solutionCheck->setEnabled(false);
//        //ui->sqrtJ2Check->setEnabled(false);
//        //ui->contourCheck->setEnabled(false);
//        //ui->polyChainCheck->setEnabled(false);
//        ui->scaleCheck->setEnabled(false);
//        ui->innerRadCheck->setEnabled(false);

//        //Disable Run buttons
//        pwbDock->ui->runBtn->setEnabled(false);
//        // disable autobar check box
//        ui->autoBar->setEnabled(false);
}

// Loads data from Wellbore analysis obj into GUI fields
void MainWindow::populate_GUIFields() {

//#ifdef PV
    TPZPlasticStepPV<TPZSandlerExtended,TPZElasticResponse> &SD = fwellb.GetCurrentConfig()->fSDPV;
    parametersDock->ui->fE->setText(QVariant(SD.fYC.E()).toString());
    parametersDock->ui->fPoisson->setText(QVariant(SD.fYC.Poisson()).toString());
    parametersDock->ui->fA->setText(QVariant(SD.fYC.A()).toString());
    parametersDock->ui->fB->setText(QVariant(SD.fYC.B()).toString());
    parametersDock->ui->fC->setText(QVariant(SD.fYC.C()).toString());
    parametersDock->ui->fD->setText(QVariant(SD.fYC.D()).toString());
    parametersDock->ui->fR->setText(QVariant(SD.fYC.R()).toString());
    parametersDock->ui->fW->setText(QVariant(SD.fYC.W()).toString());
//#else
//    TPZSandlerDimaggio<SANDLERDIMAGGIOSTEP2> &SD = fwellb.GetCurrentConfig()->fSD;
//#endif

    TPZPlasticStepPV<TPZYCMohrCoulombPV,TPZElasticResponse> &MC = fwellb.GetCurrentConfig()->fMCPV;
    parametersDock->ui->fAngleMohr->setText(QVariant(MC.fYC.Phi()).toString());
    parametersDock->ui->fCohesionMohr->setText(QVariant(MC.fYC.Cohesion()).toString());
    parametersDock->ui->fYoungMohr->setText(QVariant(MC.fYC.E()).toString());
    parametersDock->ui->fPoissonMohr->setText(QVariant(MC.fYC.Poisson()).toString());

    parametersDock->selectMat(fwellb.GetCurrentConfig()->fModel);
#ifdef DEBUG
    qDebug() << "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA fReservoirPressure : " << fwellb.GetCurrentConfig()->fReservoirPressure;
    qDebug() << "BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB fWellborePressure  : " << fwellb.GetCurrentConfig()->fWellborePressure;
#endif

    loadparamsDock->ui->Sh->setText(QVariant(fwellb.GetCurrentConfig()->fConfinementTotal.XX()).toString());
    loadparamsDock->ui->SH->setText(QVariant(fwellb.GetCurrentConfig()->fConfinementTotal.YY()).toString());
    loadparamsDock->ui->SV->setText(QVariant(fwellb.GetCurrentConfig()->fConfinementTotal.ZZ()).toString());
    pwbDock->ui->pwb->setText(QVariant(fwellb.GetCurrentConfig()->fWellborePressure).toString());
    pwbDock->ui->pres->setText(QVariant(fwellb.GetCurrentConfig()->fReservoirPressure).toString());

    REAL inner = fwellb.GetCurrentConfig()->fInnerRadius;
    REAL outer = fwellb.GetCurrentConfig()->fOuterRadius;
    int n_rad  = fwellb.GetCurrentConfig()->fNx[0]; //divisoes em X
    int n_circ = fwellb.GetCurrentConfig()->fNx[1]; //divisoes em Y
    REAL Delx = fwellb.GetCurrentConfig()->fDelx;
    geometryparamsDock->ui->r_int->setText(QVariant(inner).toString());
    geometryparamsDock->ui->r_ext->setText(QVariant(outer).toString());
    meshDock->ui->n_radial ->setText(QVariant(n_rad).toString());
    meshDock->ui->n_circ->setText(QVariant(n_circ).toString());

    prefineDock->ui->p_order->setText(QVariant(fwellb.GetCurrentConfig()->fCMesh.GetDefaultOrder()).toString());
    prefineDock->ui->isoValue->setText("0");

    REAL first_el_size = Delx/(inner*M_PI_2/n_circ);
    meshDock->ui->first_el_size->setText(QVariant(first_el_size).toString());

    //pwbDock->ui->FluidMODEL
    //Biot Coeff
    parametersDock->ui->alpha->setText(QVariant(fwellb.GetCurrentConfig()->fBiotCoef).toString());

    geometryparamsDock->setDirection(fwellb.GetCurrentConfig()->fWellConfig);
}

// Saves data from GUI into Wellbore analysis obj
void MainWindow::save_GUIFields() {

    // Setting biot coeff
    REAL biot = parametersDock->ui->alpha->text().toDouble();
    fwellb.SetBiotCoefficient(biot);

    //Seting Total WellBore Pressure and Total Confinement Stresses
    TPZManVector<STATE,3> confinement(3,0.);
    REAL ShTot = loadparamsDock->ui->Sh->text().toDouble(); //x
    REAL SHTot = loadparamsDock->ui->SH->text().toDouble(); //y
    REAL SVTot = loadparamsDock->ui->SV->text().toDouble(); //z
    REAL wellborePressure = loadparamsDock->ui->wellborePressure->text().toDouble();
    confinement[0] = ShTot;//ShEff; //sigX;
    confinement[1] = SHTot;//SHEff; //sigY;
    confinement[2] = SVTot;//SVEff; //sigZ;
    fwellb.SetConfinementTotalStresses(confinement, wellborePressure);

    // Setting reservoir Pressure // PR0
    REAL ppo = loadparamsDock->ui->porePressure->text().toDouble();
    fwellb.SetReservoirPressure(ppo);

    //Seting inner and outer radius
    REAL r_int = geometryparamsDock->ui->r_int->text().toDouble();
    REAL r_ext = geometryparamsDock->ui->r_ext->text().toDouble();
    fwellb.SetInnerOuterRadius(r_int, r_ext);

    //Seting Material Parameters
    // Sandler
    if (parametersDock->SelectedMat() == 0) {
        REAL E = parametersDock->ui->fE->text().toDouble();
        REAL poisson = parametersDock->ui->fPoisson->text().toDouble();
        REAL A = parametersDock->ui->fA->text().toDouble();
        REAL B = parametersDock->ui->fB->text().toDouble();
        REAL C = parametersDock->ui->fC->text().toDouble();
        REAL D = parametersDock->ui->fD->text().toDouble();
        REAL R = parametersDock->ui->fR->text().toDouble();
        REAL W = parametersDock->ui->fW->text().toDouble();
        fwellb.SetSanderDiMaggioParameters(poisson, E, A, B, C, R, D, W);
    }
    // Mohr
    if (parametersDock->SelectedMat() == 1) {
        REAL young = parametersDock->ui->fYoungMohr->text().toDouble();
        REAL poisson = parametersDock->ui->fPoissonMohr->text().toDouble();
        REAL angle = parametersDock->ui->fAngleMohr->text().toDouble();
        REAL cohesion = parametersDock->ui->fCohesionMohr->text().toDouble();
        fwellb.SetMohrCoulombParameters(poisson, young, cohesion, angle, angle );
    }
    // Elastic
    if (parametersDock->SelectedMat() == 2) {

        // Using very big plastic surface
        REAL young = parametersDock->ui->fYoungLinear->text().toDouble();
        REAL poisson = parametersDock->ui->fPoissonLinear->text().toDouble();

        REAL cohesion = 1.e8; // Very very big
        REAL Phi = 1.5533430342749532; // 89 degrees
        fwellb.SetMohrCoulombParameters(poisson, young, cohesion, Phi, Phi);
        fwellb.GetCurrentConfig()->fModel = EElastic;

        //fwellb.SetLinearElasticParameters(young, poisson);
    }
    //REAL A = 1; //IMPLEMENTAR: poromechDock->getA();


    //Setting geometric mesh properties
    int n_radial = meshDock->ui->n_radial->text().toInt();
    int n_circ = meshDock->ui->n_circ->text().toInt();
    int porder = prefineDock->ui->p_order->text().toInt();
    REAL first_el_size = meshDock->ui->first_el_size->text().toDouble();
    REAL Delx = first_el_size*r_int*M_PI_2/n_circ;
    TPZManVector<int,2> numdiv(2);
    numdiv[0] = n_radial; //divisoes em X
    numdiv[1] = n_circ;   //divisoes em Y
    fwellb.SetMeshTopology(Delx, numdiv);


    //Setting WellBore Configuration (**EVerticalWell** | EHorizontalWellalongh | EHorizontalWellalongH )
    fwellb.SetWellConfig(EWellConfiguration(geometryparamsDock->getDirection()));

    //Creating Computational Mesh
    fwellb.GetCurrentConfig()->CreateMesh();

    //Creating Computational PostProcMesh
    fwellb.GetCurrentConfig()->CreateComputationalMesh(porder);

    //Setting Fluid Model (ENonPenetrating - EPenetrating)
    // Fluid model // No Penetrating
    if (pwbDock->ui->FluidModelRB_NP->isChecked()) {
        fwellb.SetFluidModel(ENonPenetrating);
    }
    // Fluid model // Penetrating
    if (pwbDock->ui->FluidModelRB_P->isChecked()) {
        fwellb.SetFluidModel(EPenetrating);
    }

}

// generate a discontinuous representation of the geometric grid
vtkSmartPointer<vtkUnstructuredGrid> MainWindow::GenerateVTKGridDisc(double setZ)
{

    fwellb.GetCurrentConfig()->CreatePostProcessingMesh();
    fLoadedConfig->CreatePostProcessingMesh();

    TPZGeoMesh *gMesh = &(fLoadedConfig->fGMesh);

    TPZPostProcAnalysis &postan = fLoadedConfig->fPostprocess;

    TPZCompMesh * cmesh = postan.Mesh();

    if(!cmesh)
    {
        DebugStop();
    }

    int mirrored = ui->showFullCheck->isChecked();

    long nNodes = cmesh->NElements()*4;
    long nEls = cmesh->NElements();

    // Gerando a GRID
    vtkSmartPointer<vtkUnstructuredGrid> aGrid = vtkSmartPointer<vtkUnstructuredGrid>::New();
    // Gerando os nos / pontos
    vtkSmartPointer<vtkPoints> allPoints = vtkSmartPointer<vtkPoints>::New();
    if(!mirrored)
    {
        allPoints->SetNumberOfPoints(nNodes);
    }
    else
    {
        allPoints->SetNumberOfPoints(4*nNodes);
    }
    for (int quadrant = 0; quadrant < 4; quadrant++)
    {
        for (long el=0; el< nEls; el++) {
            TPZGeoEl *gel = cmesh->ElementVec()[el]->Reference();
            for (int locn = 0; locn<4; locn++)
            {
                int n = gel->NodeIndex(locn);
                double x = gMesh->NodeVec()[n].Coord(0);
                double y = gMesh->NodeVec()[n].Coord(1);
                if(quadrant%2) x = -x;
                if(quadrant/2) y = -y;
                double z = gMesh->NodeVec()[n].Coord(2);
                if (setZ != 0) z = setZ;
                allPoints->InsertPoint(quadrant*nNodes+el*4+locn, x, y, z);
            }
        }
        if(!mirrored) break;
    }
    aGrid->SetPoints(allPoints);

    for (int quadrant = 0; quadrant < 4; quadrant++)
    {
        // Gerando os elementos
        for (long el=0; el< nEls; el++) {
            TPZCompEl *cel = cmesh->ElementVec()[el];
                    TPZGeoEl *gel = cel->Reference();
            int gelNnodes = gel->NNodes();
            if (gelNnodes == 4 || gelNnodes == 8) { //Gerando o Quadrado
                vtkSmartPointer<vtkQuad> aQuad = vtkSmartPointer<vtkQuad>::New();
                aQuad->GetPointIds()->SetId(0, quadrant*nNodes+4*el);
                aQuad->GetPointIds()->SetId(1, quadrant*nNodes+4*el+1);
                aQuad->GetPointIds()->SetId(2, quadrant*nNodes+4*el+2);
                aQuad->GetPointIds()->SetId(3, quadrant*nNodes+4*el+3);
                aGrid->InsertNextCell(aQuad->GetCellType(), aQuad->GetPointIds());
            }
            else
            {
                DebugStop();
            }
        }
        if(!mirrored) break;
    }

//    for (int i = 0; i < ui->postPVarCBox->count(); i++) {
//        addPostProcVar(i);
//    }

    return aGrid;
}

void MainWindow::addPostProcVar(int i) {
    int mirrored = ui->showFullCheck->isChecked();

    vtkSmartPointer<vtkDoubleArray> varVtkArray = vtkSmartPointer<vtkDoubleArray>::New();

    TPZVec<REAL> varVec;

    QString var_name (ui->postPVar2CBox->itemText(i));
    QString var_type (ui->postPVar2CBox->itemData(i).toString());

    computeVisualPostProcess(varVec, var_name, var_type);

    //varVtkArray->SetNumberOfComponents(1);
    varVtkArray->SetName (var_type.toStdString().c_str());
    int varSize = varVec.size();
    int numtuples  = varSize;

    if(mirrored) numtuples *= 4;

    varVtkArray->SetNumberOfTuples(numtuples);

    // QUADRANT #:
    //
    //    1 | 0
    //   -------
    //    3 | 2
    //

    for(int quadrant =0; quadrant<4; quadrant++) {

        for (int i=0; i< varSize; i++)
        {
            double mult = 1.;

            if ( (var_name == "DisplacementMemX") ) {

                switch (quadrant) {
                    case 0:
                    mult = 1.;
                        break;
                    case 1:
                    mult = -1.;
                        break;
                    case 2:
                    mult = 1.;
                        break;
                    case 3:
                    mult = -1.;
                        break;
                }
            }
            if ( (var_name == "DisplacementMemY") ) {


                switch (quadrant) {
                    case 0:
                    mult = -1.;
                        break;
                    case 1:
                    mult = -1.;
                        break;
                    case 2:
                    mult = 1.;
                        break;
                    case 3:
                    mult = 1.;
                        break;
                }

            }

            varVtkArray->SetValue(quadrant*varSize+i, varVec[i]*mult);
        }
        if(!mirrored) break;

    }
    //aGrid->GetPointData()->SetScalars(varVtkArray);
    aGrid->GetPointData()->AddArray(varVtkArray);
#ifdef DEBUG
    qDebug() << "Creating scalar: " << varVtkArray->GetName();
#endif
}
// generate the data of var acording to the
void MainWindow::computeVisualPostProcess(TPZVec<STATE> &vecvar, QString varName, QString varType)
{
//    // Verifying exceptions
//    // Sandler
//    if (parametersDock->SelectedMat() == 0 && varName == "YieldSurface3") {
//        varName = "YieldSurface2";
//    }

    REAL coord[4][2] = {{-1.,-1.},{1.,-1.},{1.,1.},{-1.,1.}};
    TPZCompMesh *cmesh = fLoadedConfig->fPostprocess.Mesh();
TPZMaterial *mat = cmesh->MaterialVec().begin()->second;
    int varind = mat->VariableIndex(varType.toStdString());
if(varind == -1) DebugStop();
int nel = cmesh->NElements();
    vecvar.resize(4*nel);
for(int el=0; el<nel; el++)
{
TPZCompEl *cel = cmesh->ElementVec()[el];
for(int no=0; no<4; no++)
{
TPZManVector<REAL,2> xi(2),sol(1);
for(int i=0; i<2; i++) xi[i] = coord[no][i];
cel->Solution(xi,varind,sol);
//            if (    (varName == "DisplacementMemY") ||
//                    (varName == "PrincipalStress2") ||
//                    (varName == "PrincipalStrain2") ||
//                    (varName == "YieldSurface2") ||
//                    (varName == "NormalStressY") ||
//                    (varName == "NormalStrainY") ||
//                    (varName == "ShearStressXZ") ||
//                    (varName == "ShearStrainXZ") ||
//                    (varName == "TotalPlasticStrainY"))
//                vecvar[4*el+no] = IsZero(sol[1])? 0.0 : (float)sol[1]; //removing double precision
//            else
//            if (    (varName == "DisplacementMemZ") ||
//                    (varName == "PrincipalStress3") ||
//                    (varName == "PrincipalStrain3") ||
//                    (varName == "YieldSurface3") ||
//                    (varName == "NormalStressZ") ||
//                    (varName == "NormalStrainZ") ||
//                    (varName == "ShearStressYZ") ||
//                    (varName == "ShearStrainYZ") ||
//                    (varName == "TotalPlasticStrainZ"))
//                vecvar[4*el+no] = IsZero(sol[2])? 0.0 : (float)sol[2]; //removing double precision
//            else
                //default for scalar and first position from vecs
                vecvar[4*el+no] = IsZero(sol[0])? 0.0 : (float)sol[0]; //removing double precision
}
}
}

void MainWindow::composeScene() {

    //Creating mappers
    mapperMesh = vtkSmartPointer<vtkDataSetMapper>::New();
    mapperWire = vtkSmartPointer<vtkDataSetMapper>::New();
    mapperEllips = vtkSmartPointer<vtkPolyDataMapper>::New();
    mapperPolychain = vtkSmartPointer<vtkPolyDataMapper>::New();
    contMapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    contMapper2 = vtkSmartPointer<vtkPolyDataMapper>::New();
    mapperInnerRad = vtkSmartPointer<vtkPolyDataMapper>::New();

    //Creating Grid
    aGrid = vtkSmartPointer<vtkUnstructuredGrid>::New();

    //Creating actors
    actorMesh = vtkSmartPointer<vtkActor>::New();
    actorWire = vtkSmartPointer<vtkActor>::New();
    actorEllips = vtkSmartPointer<vtkActor>::New();
    actorPolychain = vtkSmartPointer<vtkActor>::New();
    actorCont = vtkSmartPointer<vtkActor>::New();
    actorCont2 = vtkSmartPointer<vtkActor>::New();
    actorInnerRad = vtkSmartPointer<vtkActor>::New();

    //Creating isolines
    contours = vtkSmartPointer<vtkContourFilter>::New();
    contours2 = vtkSmartPointer<vtkContourFilter>::New();

    //Creating scalar bar
    scalarBar = vtkSmartPointer<vtkScalarBarActor>::New();

    //Creating scene and attaching to QVtk component
    renderer = vtkSmartPointer<vtkRenderer>::New();
    renderer->GetActiveCamera()->ParallelProjectionOn();
    renderer_window = vtkSmartPointer<vtkRenderWindow>::New();
    renderer_window->AddRenderer(renderer);
    interactor = ui->qvtkWidget->GetInteractor();
    vtkSmartPointer<WellboreInteractorStyle> style = vtkSmartPointer<WellboreInteractorStyle>::New();
    interactor->SetInteractorStyle( style );
    interactor->SetRenderWindow(renderer_window);
    ui->qvtkWidget->SetRenderWindow(renderer_window);

    //Adding actors to the scene
    //renderer->AddActor2D(scalarBar);
    renderer->AddActor(actorEllips);
    renderer->AddActor(actorPolychain);
    renderer->AddActor(actorMesh);
    renderer->AddActor(actorWire);
    renderer->AddActor(actorCont);
    renderer->AddActor(actorCont2);
    renderer->AddActor(actorInnerRad);

    initialZoomVTK();
}

void MainWindow::openBin() {
    // Gerando a GRID
    aGrid = this->GenerateVTKGridDisc(-1.);

    postProcessVarChanged();

    //Updating UI - enabling checkboxes
    ui->gridCheck->setEnabled(true);
    ui->solutionCheck->setEnabled(true);
    ui->innerRadCheck->setEnabled(true);

    ui->showFullCheck->setEnabled(true);
    ui->scaleCheck->setEnabled(true);

    ui->autoBar->setEnabled(true);

    changedTabs_updateCheckboxes(true);	
}

void MainWindow::initialZoomVTK() {
    // ZOOM fixes
    renderer->ResetCamera();
    renderer_window->Render();

    //renderer->GetActiveCamera()->SetParallelProjection(true);
    REAL inner = fLoadedConfig->fInnerRadius;
    //renderer->GetActiveCamera()->SetPosition(0.7*inner,0,1);
    //renderer->GetActiveCamera()->SetFocalPoint(0.7*inner,0,0);
    double *position = renderer->GetActiveCamera()->GetPosition();
    renderer->GetActiveCamera()->SetPosition(inner*1.5,inner/2.2,position[2]);
    renderer->GetActiveCamera()->SetFocalPoint(inner*1.5,inner/2.2,0.);
    renderer->GetActiveCamera()->SetParallelScale(inner/2.);
    //vtkIndent indent;
    //renderer->PrintSelf(std::cout,indent);
    //renderer->GetActiveCamera()->PrintSelf(std::cout,indent);
    renderer_window->Render();
}

void MainWindow::BreakApplyBtn_clicked()
{
    int configCount = fwellb.GetConfigListSize();

    if (configCount > 0) {

        int logrow = writeHistory2(tr("Applying breakout..."));

        // Update a & b values
        REAL a = identifyellipsDock->ui->fa->text().toDouble();
        REAL b = identifyellipsDock->ui->fb->text().toDouble();

        REAL sqrtJ2 = identifyellipsDock->ui->sqrtJ2->text().toDouble();

        fwellb.AddEllipticBreakout(a, b);
        fwellb.GetCurrentConfig()->ModifyWellElementsToQuadratic();
        fwellb.GetCurrentConfig()->ComputeElementDeformation();
        fwellb.GetCurrentConfig()->CreatePostProcessingMesh();

        //Saving TConfig
        std::stringstream strout;
        strout << "Breakout sqrtJ2 = " << sqrtJ2;
        fwellb.SaveConfig(strout);

        reloadHistory();

        writeHistory2(tr("...done"), true, logrow);

        fLoadedConfig = fwellb.GetCurrentConfig();

        openBin();
    }
#ifdef DEBUG
    qDebug() << "on_BreakApplyBtn_clicked ok";
#endif
}

void MainWindow::divideBtn_clicked()
{
    double isovalue = hrefineDock->ui->isoValue->text().toDouble();

    int logrow = writeHistory2(QString(tr("Dividing elements above ")).append(QVariant(isovalue).toString()));

    unsigned int numelements;
    numelements = fwellb.DivideElementsAbove(isovalue);

    //Saving TConfig
    std::stringstream strout;
    strout << "Divide H sqrtJ2 <= " << isovalue ;
    fwellb.SaveConfig(strout);

    QString out = QString(tr(" Number of elements divided ")).append(QVariant(numelements).toString());

    writeHistory2(out, true, logrow);

    reloadHistory();

    fLoadedConfig = fwellb.GetCurrentConfig();
    
    openBin();

#ifdef DEBUG
    qDebug() << "on_divideBtn_clicked ok";
#endif
}

void MainWindow::prefineBtn_clicked()
{
    double isovalue = prefineDock->ui->isoValue->text().toDouble();
    int p_order = prefineDock->ui->p_order->text().toInt();

    int logrow = writeHistory2(tr("P refining..."));

    unsigned int numelements;
    numelements = fwellb.PRefineElementAbove(isovalue, p_order);

    //Saving TConfig
    std::stringstream strout;
    strout << "Refine P=" << p_order << " sqrtJ2 <= " << isovalue;
    fwellb.SaveConfig(strout);

    QString out = QString(tr(" Number of elements prefined ")).append(QVariant(numelements).toString());
    writeHistory2(out, true, logrow);

    reloadHistory();

    fLoadedConfig = fwellb.GetCurrentConfig();

    openBin();

#ifdef DEBUG
    qDebug() << "on_prefineBtn_clicked ok";
#endif
}

void MainWindow::calculateABBtn_clicked()
{
    int configCount = fwellb.GetConfigListSize();

    if (configCount > 0) {

        int logrow = writeHistory2(tr("Calculating A and B..."));

        // Update contour line for SQRT J2 value from Interface
        REAL value_sqrtJ2 = identifyellipsDock->ui->sqrtJ2->text().toDouble();

        REAL a,b;

        fwellb.ComputeAandB(value_sqrtJ2, a,b);

        // Update a & b values
        //
        identifyellipsDock->ui->fa->setText(QString("%1").arg(a, 0, 'g', 5));
        identifyellipsDock->ui->fb->setText(QString("%1").arg(b, 0, 'g', 5));

        // Update ellipse curve on VTK screen
        redrawEllips();

        redrawPolyChain();

        writeHistory2(tr("...done"), true, logrow);
    }

#ifdef DEBUG
    qDebug() << "on_calculateABBtn_clicked ok";
#endif
}

void MainWindow::redrawJ2() {
    // Update contour line for SQRT J2 value from Interface
    REAL value_sqrtJ2 = identifyellipsDock->ui->sqrtJ2->text().toDouble();
#ifdef DEBUG
    std::cout << __PRETTY_FUNCTION__ << "value sqrtJ2 " << value_sqrtJ2 << std::endl;
#endif
#if VTK_MAJOR_VERSION <= 5
    contours2->SetInput(aGrid);
#else
    contours2->SetInputData(aGrid);
#endif
    contours2->GenerateValues(1, value_sqrtJ2, value_sqrtJ2);//(1, 0.0001, 0.0001);
    //contours2->ComputeScalarsOff();
    // map the contours to graphical primitives
#if VTK_MAJOR_VERSION <= 5
    contMapper2->SetInput(contours2->GetOutput());
#else
    contMapper2->SetInputData(contours2->GetOutput());
#endif
    //contMapper2->SetScalarRange(0.0001, 0.0001);
    contMapper2->ScalarVisibilityOff();
    // Update an actor for the contours
    actorCont2->SetMapper(contMapper2);
    actorCont2->GetProperty()->SetColor(0., 0., 0.);
    actorCont2->GetProperty()->SetLineWidth(3);

    renderer_window->Render();
}

void MainWindow::redrawPolyChain() {
    actorPolychain->SetMapper(mapperPolychain);
    actorPolychain->GetProperty()->SetColor(0., 0., 1.);
    actorPolychain->GetProperty()->SetLineWidth(2.5);
    actorPolychain->GetProperty()->SetPointSize(5);

    vtkSmartPointer<vtkPolyData> pd = vtkSmartPointer<vtkPolyData>::New();

    vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
    vtkSmartPointer<vtkCellArray> verts = vtkSmartPointer<vtkCellArray>::New();

    int mirrored = ui->showFullCheck->isChecked();

    // Update polychain line based on SQRT J2 value from Interface
    REAL value_sqrtJ2 = identifyellipsDock->ui->sqrtJ2->text().toDouble();

#ifdef DEBUG
    std::cout << __PRETTY_FUNCTION__ << "value sqrtJ2 " << value_sqrtJ2 << std::endl;
#endif

    std::multimap<REAL, REAL> polygonalChainbase;
    fwellb.GetJ2Isoline(value_sqrtJ2, polygonalChainbase);
    //REAL maxy = fwellb.GetCurrentConfig()->MaxYfromLastBreakout();
    REAL maxy = fLoadedConfig->MaxYfromLastBreakout();

    vtkSmartPointer<vtkIdList> lineIndices = vtkSmartPointer<vtkIdList>::New();

    // Setup colors
    vtkSmartPointer<vtkUnsignedCharArray> colors = vtkSmartPointer<vtkUnsignedCharArray>::New();
    colors->SetNumberOfComponents(3);
    colors->SetName ("Colors");
    unsigned char red[3] = {255, 0, 0};
    unsigned char blue[3] = {0, 0, 255};

    int i=0;
    for (int quadrant = 0; quadrant < 4; quadrant++)
    {
        for (std::multimap<REAL, REAL>::iterator it = polygonalChainbase.begin(); it != polygonalChainbase.end(); it++) {
            double x = it->first;
            double y = it->second;
            if(quadrant%2) x = -x;
            if(quadrant/2) y = -y;

            points->InsertPoint(static_cast<vtkIdType>(i), x, y, 0.0 );
            lineIndices->InsertNextId(static_cast<vtkIdType>(i));

            TPZManVector<REAL,2> co(2);
            co[0] = 0.;
            co[1] = it->second;
            //fwellb.GetCurrentConfig()->ProjectNode(co);
            fLoadedConfig->ProjectNode(co);

            //Color
            //if (it->second < maxy && it->first-co[0] > 0.01*(fwellb.GetCurrentConfig()->fInnerRadius))
            if (it->second < maxy && it->first-co[0] > 0.01*(fLoadedConfig->fInnerRadius))
                colors->InsertNextTupleValue (blue);
            else
                colors->InsertNextTupleValue (red);

            i++;
        }
        if(!mirrored) break;
    }

    verts->InsertNextCell(lineIndices);

    pd->SetPoints(points);
    pd->SetVerts(verts);
    pd->Squeeze();

    // Setup colors
    pd->GetPointData()->SetScalars(colors);

    // UPDATE the polychain curve with new points
#if VTK_MAJOR_VERSION <= 5
    mapperPolychain->SetInput(pd);
#else
    mapperPolychain->SetInputData(pd);
#endif

    renderer_window->Render();
}

void MainWindow::redrawEllips() {
    actorEllips->SetMapper(mapperEllips);
    actorEllips->GetProperty()->SetColor(0., 1., 0.); //green
    actorEllips->GetProperty()->SetLineWidth(2.5);

    vtkSmartPointer<vtkPolyData> pd = vtkSmartPointer<vtkPolyData>::New();

    vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
    vtkSmartPointer<vtkCellArray> lines = vtkSmartPointer<vtkCellArray>::New();

    int mirrored = ui->showFullCheck->isChecked();

    int CenterX , CenterY;
    CenterX = 0;
    CenterY = 0;
    // Update a & b values
    REAL r1 = 0;
    REAL r2 = 0;
    REAL sh = loadparamsDock->ui->Sh->text().toDouble();
    REAL SH = loadparamsDock->ui->SH->text().toDouble();
    if (sh > SH) {
        r1 = identifyellipsDock->ui->fa->text().toDouble();
        r2 = identifyellipsDock->ui->fb->text().toDouble();
    }
    else
    {
        r2 = identifyellipsDock->ui->fa->text().toDouble();
        r1 = identifyellipsDock->ui->fb->text().toDouble();
    }

    vtkSmartPointer<vtkIdList> lineIndices = vtkSmartPointer<vtkIdList>::New();

    // number of points in ellipse curve
    const int pontos=50;
    for (int i = 0; i< pontos; i++)
    {
        double quad = 2.0*vtkMath::Pi();
        if (!mirrored) quad /= 4;
        double angle = quad*i/(float(pontos)-1);
        points->InsertPoint(static_cast<vtkIdType>(i), r1 * cos(angle) + CenterX, r2 * sin(angle) + CenterY, 0.0 );
        lineIndices->InsertNextId(static_cast<vtkIdType>(i));
    }

    lines->InsertNextCell(lineIndices);

    pd->SetPoints(points);
    pd->SetLines(lines);
    pd->Squeeze();

    // UPDATE the ellips curve with new points
#if VTK_MAJOR_VERSION <= 5
    mapperEllips->SetInput(pd);
#else
    mapperEllips->SetInputData(pd);
#endif

    renderer_window->Render();
}

void MainWindow::redrawInnerRad() {
    actorInnerRad->SetMapper(mapperInnerRad);
    actorInnerRad->GetProperty()->SetColor(1., 1., 1.); //white
    actorInnerRad->GetProperty()->SetLineWidth(2.5);

    vtkSmartPointer<vtkPolyData> pd = vtkSmartPointer<vtkPolyData>::New();

    vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
    vtkSmartPointer<vtkCellArray> lines = vtkSmartPointer<vtkCellArray>::New();

    int mirrored = ui->showFullCheck->isChecked();

    // Get inner rad value
    REAL innerRad = geometryparamsDock->ui->r_int->text().toDouble();

    vtkSmartPointer<vtkIdList> lineIndices = vtkSmartPointer<vtkIdList>::New();

    // number of points in radius curve
    const int pontos=50;
    for (int i = 0; i< pontos; i++)
    {
        double quad = 2.0*vtkMath::Pi();
        if (!mirrored) quad /= 4;
        double angle = quad*i/(float(pontos)-1);
        //points->InsertPoint(static_cast<vtkIdType>(i), r1 * cos(angle) + CenterX, r2 * sin(angle) + CenterY, 0.0 );
        points->InsertPoint(static_cast<vtkIdType>(i), innerRad * cos(angle), innerRad * sin(angle), 0.0 );
        lineIndices->InsertNextId(static_cast<vtkIdType>(i));
    }

    lines->InsertNextCell(lineIndices);

    pd->SetPoints(points);
    pd->SetLines(lines);
    pd->Squeeze();

    // UPDATE the radius curve with new points
#if VTK_MAJOR_VERSION <= 5
    mapperInnerRad->SetInput(pd);
#else
    mapperInnerRad->SetInputData(pd);
#endif

    renderer_window->Render();
}

void MainWindow::on_gridCheck_toggled(bool checked)
{
#ifdef DEBUG
    qDebug() << "on_gridCheck_toggled " << checked;
#endif

    if (actorWire.GetPointer() != NULL) {
        actorWire->SetVisibility(checked);
        renderer_window->Render();
    }
}

void MainWindow::on_solutionCheck_toggled(bool checked)
{
#ifdef DEBUG
    qDebug() << "on_solutionCheck_toggled " << checked;
#endif

    if (actorMesh.GetPointer() != NULL) {
        actorMesh->SetVisibility(checked);
        renderer_window->Render();
    }
}

void MainWindow::ellipsCheck_toggled(bool checked)
{
#ifdef DEBUG
    qDebug() << "on_ellipsCheck_toggled " << checked;
#endif

    if (actorEllips.GetPointer() != NULL) {
        actorEllips->SetVisibility(checked);
        renderer_window->Render();
    }
}

void MainWindow::sqrtJ2Check_toggled(bool checked)
{
#ifdef DEBUG
    qDebug() << "on_sqrtJ2Check_toggled " << checked;
#endif

    if (actorCont2.GetPointer() != NULL) {
        actorCont2->SetVisibility(checked);
        renderer_window->Render();
    }
}

void MainWindow::contourCheck_toggled(bool checked)
{
#ifdef DEBUG
    qDebug() << "on_contourCheck_toggled " << checked;
#endif

    if (actorCont.GetPointer() != NULL) {
        actorCont->SetVisibility(checked);
        renderer_window->Render();
    }
}

void MainWindow::on_showFullCheck_toggled(bool checked)
{
#ifdef DEBUG
    qDebug() << "on_showFullCheck_toggled " << checked;
#endif
    openBin();
}

void MainWindow::polyChainCheck_toggled(bool checked)
{
#ifdef DEBUG
    qDebug() << "on_polyChainCheck_toggled " << checked;
#endif

    if (actorPolychain.GetPointer() != NULL) {
        actorPolychain->SetVisibility(checked);
        renderer_window->Render();
    }
}

void MainWindow::on_innerRadCheck_toggled(bool checked)
{
#ifdef DEBUG
    qDebug() << "on_innerRadCheck_toggled " << checked;
#endif

    if (actorInnerRad.GetPointer() != NULL) {
        actorInnerRad->SetVisibility(checked);
        renderer_window->Render();
    }
}

void MainWindow::on_resetZoomBtn_clicked()
{
    initialZoomVTK();
}

void MainWindow::on_actionInitialConfig_triggered(bool checked)
{
    initialDialog->show();
}

void MainWindow::on_scaleCheck_toggled(bool checked)
{
#ifdef DEBUG
    qDebug() << "on_scaleCheck_toggled " << checked;
#endif

    if (scalarBar.GetPointer() != NULL) {
        scalarBar->SetVisibility(checked);
        renderer_window->Render();
    }
}

void MainWindow::pwb_RunBtn() {

    pwbDock->ui->runBtn->setEnabled(false);

    QString file_basename ("WellboreQT");
    QString file_extname (".bin");
    std::cout << std::setprecision(15);

    int configCount = fwellb.GetConfigListSize();

    int logrow2 = writeHistory2(QString(tr("Running step ")).append(QVariant(configCount).toString()).append("..."));

    if (configCount == 0) {

        //Creating Post-ProcessMesh
        fwellb.GetCurrentConfig()->CreatePostProcessingMesh();

        QString filename (file_basename);
        filename.append(QVariant(configCount).toString()).append(file_extname);

    #ifdef TOBIN
        TPZBFileStream save;
        save.OpenWrite(filename.toStdString());
        fwellb.Write(save);
    #endif

        std::stringstream strout;
        strout << "Initial Config";
        fwellb.SaveConfig(strout);

        reloadHistory();

    }
    else if (configCount == 1) {
        //Setting Fluid Model (ENonPenetrating - EPenetrating)
        // Fluid model // No Penetrating
        if (pwbDock->ui->FluidModelRB_NP->isChecked()) {
            fwellb.SetFluidModel(ENonPenetrating);
        }
        // Fluid model // Penetrating
        if (pwbDock->ui->FluidModelRB_P->isChecked()) {
            fwellb.SetFluidModel(EPenetrating);
        }

        int num_steps = pwbDock->ui->steps->text().toInt();
        int numnewton = 80;

        fwellb.GetCurrentConfig()->ModifyWellElementsToQuadratic();
        fwellb.ExecuteInitialSimulation(num_steps, numnewton);
        reloadHistory();
    }
    else
    {
        // Fluid model // No Penetrating
        if (pwbDock->ui->FluidModelRB_NP->isChecked()) {
            fwellb.SetFluidModel(ENonPenetrating);
        }
        // Fluid model // Penetrating
        if (pwbDock->ui->FluidModelRB_P->isChecked()) {
            fwellb.SetFluidModel(EPenetrating);
        }
        
        int nsteps = pwbDock->ui->steps->text().toInt();
        REAL pwb = pwbDock->ui->pwb->text().toDouble();
        //REAL pressure = fwellb.GetCurrentConfig()->fFluidPressure; // current pressure
        REAL pressure = fwellb.GetCurrentConfig()->fWellborePressure; // current pressure
        REAL increment = (pressure - pwb) / nsteps;
        
#ifdef DEBUG
        qDebug() << " #################################" << pressure << "+++++++++++++" << increment;
#endif
        REAL ppo = pwbDock->ui->pres->text().toDouble();
        
        fwellb.EvolveBothPressures(nsteps, pwb, ppo);
        
        QString filename (file_basename);
        filename.append(QVariant(configCount).toString()).append(file_extname);
        reloadHistory();
#ifdef TOBIN
        TPZBFileStream save;
        save.OpenWrite(filename.toStdString());
        fwellb.Write(save);
#endif
        
        reloadHistory();

    }

    openBin();

    writeHistory2(tr("...done"), true, logrow2);

#ifdef DEBUG
    qDebug() << "pwb_RunBtn() ok";
#endif
    pwbDock->ui->runBtn->setEnabled(true);
}

void MainWindow::historyClicked(int row) {
#ifdef DEBUG
    qDebug() << "historyClicked row=" << row ;
#endif
    if (row == -1) {
#ifdef DEBUG
        qDebug() << " ERRO historyClicked ";
#endif
        //DebugStop();
        return;
    }

    // GET CONFIG INDEX STORED AT DATA-ROLE=5 FROM HISTORY ITEM
    row = historyDock->ui->historyList->item(row)->data(5).toInt();

#ifdef DEBUG
    qDebug() << "historyClicked config idx=" << row ;
#endif

    if (row >= fwellb.GetConfigListSize()) {
        fLoadedConfig = fwellb.GetCurrentConfig();
#ifdef DEBUG
        qDebug() << "CURRENT ***********************************************";
#endif
    }
    else
        fLoadedConfig = fwellb.GetConfig(row);

    //loadData();
    openBin();

#ifdef DEBUG
    qDebug() << "historyClicked ok";
#endif
}

void MainWindow::changedMinMax() {
    // save values
    minVal = ui->minBar->text().toDouble();
    maxVal = ui->maxBar->text().toDouble();

    // redraw using new values
    postProcessVarChanged ();
}

void MainWindow::postProcessVarChanged () {

    if (ui->postPVar2CBox->currentIndex() == -1)
        return;
    if (ui->postPVar1CBox->currentIndex() == -1)
        return;

#ifdef DEBUG
    qDebug() << "NEWpostProcessVarChanged " << " current1=" <<ui->postPVar1CBox->currentIndex()
                                            << " current2=" <<ui->postPVar2CBox->currentIndex()
                                            << " value=" <<ui->postPVar2CBox->itemData(ui->postPVar2CBox->currentIndex()).toString();
#endif

    addPostProcVar(ui->postPVar2CBox->currentIndex());

    aGrid->GetPointData()->SetActiveScalars(ui->postPVar2CBox->itemData(ui->postPVar2CBox->currentIndex()).toString().toStdString().c_str());

    // Gerando os nos / pontos
#if VTK_MAJOR_VERSION <= 5
    mapperMesh->SetInput(aGrid);
    mapperWire->SetInput(aGrid);
#else
    mapperMesh->SetInputData(aGrid);
    mapperWire->SetInputData(aGrid);
#endif
    mapperWire->ScalarVisibilityOff();

    if (ui->autoBar->isChecked()) {
        double *scalar_range = aGrid->GetScalarRange();
        //changed displayed min / max values, BUT DONT SAVE THEM
        ui->minBar->setText(QString("%1").arg(scalar_range[0], 0, 'g', 5));
        ui->maxBar->setText(QString("%1").arg(scalar_range[1], 0, 'g', 5));
#ifdef DEBUG
        qDebug() << "AUTO_RANGE: " << scalar_range[0] << " ____ " << scalar_range[1];
        qDebug() << "SAVED_RANGE: " << minVal << " ____ " << maxVal;
#endif
        mapperMesh->SetScalarRange(scalar_range);
    }
    else
    {
        // load min / max saved values
        double scalar_range[2];
        scalar_range[0] = minVal;
        scalar_range[1] = maxVal;
        ui->minBar->setText(QString("%1").arg(scalar_range[0], 0, 'g', 5));
        ui->maxBar->setText(QString("%1").arg(scalar_range[1], 0, 'g', 5));
#ifdef DEBUG
        qDebug() << "ENTERED_RANGE: " << scalar_range[0] << " ____ " << scalar_range[1];
        qDebug() << "SAVED_RANGE: " << minVal << " ____ " << maxVal;
#endif
        mapperMesh->SetScalarRange(scalar_range);
    }


    actorMesh->SetMapper(mapperMesh);
    actorMesh->PickableOff();
    //actorMesh->GetProperty()->SetOpacity(0.5);

    actorWire->SetMapper(mapperWire);
    actorWire->GetProperty()->SetRepresentationToWireframe();
    actorWire->GetProperty()->SetColor(0., 0., 0.);

    vtkSmartPointer<vtkLookupTable> lut = vtkSmartPointer<vtkLookupTable>::New();
//    int tableSize; = 256; //std::max(resolution*resolution + 1, 10);
//    lut->SetNumberOfTableValues(tableSize);
    //lut->SetRange(1, 10);
    //lut->SetValueRange(1, 5);
//    lut->SetNumberOfTableValues(1);
//    lut->set
//    lut->SetTableRange (0, 1);
//    lut->SetHueRange (0, 1);
//    lut->SetSaturationRange (1, 1);
//    lut->SetValueRange (1, 1);
    lut->Build();

    // Updating scalarBar
    // Adding scalarBar actor if it is not there yet:
    //     this fixes the Warning: vtkScalarBarActor: Need a mapper to render a scalar bar
    if (! renderer->GetActors2D()->IsItemPresent(scalarBar))
        renderer->AddActor2D(scalarBar);
    ////////////////////////scalarBar->SetLookupTable(mapperMesh->GetLookupTable());
    mapperMesh->SetLookupTable(lut);
    scalarBar->SetLookupTable(lut);
    scalarBar->SetTitle(ui->postPVar2CBox->itemData(ui->postPVar2CBox->currentIndex()).toString().toStdString().c_str());
    vtkSmartPointer<vtkTextProperty> titleProp = vtkSmartPointer<vtkTextProperty>::New();
    titleProp->SetColor(1.,1.,1.);
    titleProp->SetFontFamily(VTK_ARIAL);
    titleProp->SetFontSize(10);
    titleProp->SetBold(1);
    scalarBar->SetTitleTextProperty(titleProp);
    vtkSmartPointer<vtkTextProperty> labelProp = vtkSmartPointer<vtkTextProperty>::New();
    labelProp->SetColor(1.,1.,1.); //scalarBar->GetLabelTextProperty()->SetColor(1., 1., 1.);
    labelProp->SetFontFamily(VTK_ARIAL);
    labelProp->SetFontSize(10);
    labelProp->SetBold(1);
    scalarBar->SetLabelTextProperty(labelProp);
    //scalarBar->SetNumberOfLabels(4);
    scalarBar->GetPositionCoordinate()->SetCoordinateSystemToNormalizedViewport();
    scalarBar->GetPositionCoordinate()->SetValue(0.01, 0.1);
    scalarBar->SetOrientationToVertical();
    scalarBar->SetWidth(0.1);
    scalarBar->SetHeight(0.7);
    scalarBar->SetNumberOfLabels(4);

    // Update contour lines based on stored solutions
#if VTK_MAJOR_VERSION <= 5
    contours->SetInput(aGrid);
#else
    contours->SetInputData(aGrid);
#endif
    contours->GenerateValues(10, 0.0001, 0.007);
    //contours->ComputeScalarsOff();
    // map the contours to graphical primitives
#if VTK_MAJOR_VERSION <= 5
    contMapper->SetInput(contours->GetOutput());
#else
    contMapper->SetInputData(contours->GetOutput());
#endif
    //contMapper->SetScalarRange(0.0001, 0.001);
    contMapper->ScalarVisibilityOff();
    // Update an actor for the contours
    actorCont->SetMapper(contMapper);
    actorCont->GetProperty()->SetColor(1., 1., 1.);
    actorCont->GetProperty()->SetLineWidth(1);

    // SQRT J2
    // Update contour line for SQRT J2 value from Interface
    redrawJ2();

    // Recalculate and redraw ellipse
    redrawEllips();

    // Redraw Polychain
    redrawPolyChain();

    //Redraw inner radius
    redrawInnerRad();

    renderer_window->Render();

#ifdef DEBUG
    qDebug() << "postProcessVarChanged ok";
#endif
}

void MainWindow::on_actionRestart_triggered(bool checked)
{
#ifdef DEBUG
    qDebug() << "on_actionRestart_triggered " << checked;
#endif

    TPZWellBoreAnalysis a;

    fwellb = a;


    fLoadedConfig = fwellb.GetCurrentConfig();

    loadDefaultData();

    populate_GUIFields();

    composeScene();

#ifdef DEBUG
    qDebug() << "on_actionRestart_triggered - ok ";
#endif
}

void MainWindow::fluidModelChanged (int fluidModel) {
#ifdef DEBUG
    qDebug() << "fluidModelChanged " << fluidModel;
#endif

    fwellb.SetFluidModel(EFluidModel(fluidModel));
}

void MainWindow::on_actionSave_triggered()
{
#ifdef DEBUG
    qDebug() << "saveSimulation";
#endif

    QString fileName = QFileDialog::getSaveFileName(this,
        tr("Save Simulation"), "",
        tr("SEST2D simulation files (*.ses);;All Files (*)"));

    if (fileName.isEmpty())
        return;
    else {
        TPZBFileStream save;
        save.OpenWrite(fileName.toStdString());
        fwellb.Write(save);
    }
}


//        QFile file(fileName);

//        if (!file.open(QIODevice::ReadOnly)) {
//            QMessageBox::information(this, tr("Unable to open file"),
//                file.errorString());
//            return;
//        }

//        QDataStream in(&file);
//        in.setVersion(QDataStream::Qt_4_5);
//        contacts.empty();   // empty existing contacts
//        in >> contacts;

//        if (contacts.isEmpty()) {
//            QMessageBox::information(this, tr("No contacts in file"),
//                tr("The file you are attempting to open contains no contacts."));
//        } else {
//            QMap<QString, QString>::iterator i = contacts.begin();
//            nameLine->setText(i.key());
//            addressText->setText(i.value());
//        }

void MainWindow::on_actionOpen_triggered()
{
#ifdef DEBUG
    qDebug() << "loadSimulation";
#endif

    QString fileName = QFileDialog::getOpenFileName(this,
        tr("Open Simulation"), "",
        tr("SEST2D simulation files (*.ses);;All Files (*)"));

    if (fileName.isEmpty())
        return;
    else {

        TPZWellBoreAnalysis a;
        fwellb = a;

        TPZBFileStream read;
        read.OpenRead(fileName.toStdString());
        fwellb.Read(read);


        fLoadedConfig = fwellb.GetCurrentConfig();

        composeScene();

        openBin();

        populate_GUIFields();

        reloadHistory();

        initialZoomVTK();

    }
}

void MainWindow::populate_postPVar2CBox (int idx) {
    //qDebug() << "Entrou populate_postPVar2CBox";

    int c1_idx = ui->postPVar1CBox->currentIndex();

    ui->postPVar2CBox->clear();

    switch (c1_idx) {
    case 0:
        ui->postPVar2CBox->addItem("Vol","StrainVol");
        ui->postPVar2CBox->addItem("XX","StrainXX");
        ui->postPVar2CBox->addItem("YY","StrainYY");
        ui->postPVar2CBox->addItem("ZZ","StrainZZ");
        ui->postPVar2CBox->addItem("XY","StrainXY");
        ui->postPVar2CBox->addItem("XZ","StrainXZ");
        ui->postPVar2CBox->addItem("YZ","StrainYZ");
    break;
    case 1:
        ui->postPVar2CBox->addItem("Vol","ElStrainVol");
        ui->postPVar2CBox->addItem("XX","ElStrainXX");
        ui->postPVar2CBox->addItem("YY","ElStrainYY");
        ui->postPVar2CBox->addItem("ZZ","ElStrainZZ");
        ui->postPVar2CBox->addItem("XY","ElStrainXY");
        ui->postPVar2CBox->addItem("XZ","ElStrainXZ");
        ui->postPVar2CBox->addItem("YZ","ElStrainYZ");
    break;
    case 2:
        ui->postPVar2CBox->addItem("Vol","PlStrainVol");
        ui->postPVar2CBox->addItem("XX","PlStrainXX");
        ui->postPVar2CBox->addItem("YY","PlStrainYY");
        ui->postPVar2CBox->addItem("ZZ","PlStrainZZ");
        ui->postPVar2CBox->addItem("XY","PlStrainXY");
        ui->postPVar2CBox->addItem("XZ","PlStrainXZ");
        ui->postPVar2CBox->addItem("YZ","PlStrainYZ");
        ui->postPVar2CBox->addItem("PlasticSqJ2","PlStrainSqJ2");
        ui->postPVar2CBox->addItem("PlasticSqJ2El","PlStrainSqJ2El");
        ui->postPVar2CBox->addItem("Alpha","PlAlpha");
    break;
    case 3:
        ui->postPVar2CBox->addItem("X","DisplacementX");
        ui->postPVar2CBox->addItem("Y","DisplacementY");
        ui->postPVar2CBox->addItem("Z","DisplacementZ");
        ui->postPVar2CBox->addItem("Total","DisplacementTotal");
    break;
    case 4:
        ui->postPVar2CBox->addItem("POrder","POrder");
        ui->postPVar2CBox->addItem("Steps","NSteps");
    break;
    case 5:
        ui->postPVar2CBox->addItem("I1","TotStressI1");
        ui->postPVar2CBox->addItem("J2","TotStressJ2");
        ui->postPVar2CBox->addItem("XX","TotStressXX");
        ui->postPVar2CBox->addItem("YY","TotStressYY");
        ui->postPVar2CBox->addItem("ZZ","TotStressZZ");
        ui->postPVar2CBox->addItem("XY","TotStressXY");
        ui->postPVar2CBox->addItem("XZ","TotStressXZ");
        ui->postPVar2CBox->addItem("YZ","TotStressYZ");
        ui->postPVar2CBox->addItem("1","TotStress1");
        ui->postPVar2CBox->addItem("2","TotStress2");
        ui->postPVar2CBox->addItem("3","TotStress3");
        ui->postPVar2CBox->addItem("ETotStressRR");
        ui->postPVar2CBox->addItem("ETotStressRT");
        ui->postPVar2CBox->addItem("ETotStressTT");
    break;
    case 6:
        ui->postPVar2CBox->addItem("Pore Pressure","PorePressure");
    break;
    case 7:
        ui->postPVar2CBox->addItem("I1","EffStressI1");
        ui->postPVar2CBox->addItem("J2","EffStressJ2");
        ui->postPVar2CBox->addItem("XX","EffStressXX");
        ui->postPVar2CBox->addItem("YY","EffStressYY");
        ui->postPVar2CBox->addItem("ZZ","EffStressZZ");
        ui->postPVar2CBox->addItem("XY","EffStressXY");
        ui->postPVar2CBox->addItem("XZ","EffStressXZ");
        ui->postPVar2CBox->addItem("YZ","EffStressYZ");
        ui->postPVar2CBox->addItem("1","EffStress1");
        ui->postPVar2CBox->addItem("2","EffStress2");
        ui->postPVar2CBox->addItem("3","EffStress3");
    break;
    case 8:
        ui->postPVar2CBox->addItem("1","YieldSurface1");
        ui->postPVar2CBox->addItem("2","YieldSurface2");
        ui->postPVar2CBox->addItem("3","YieldSurface3");
    break;
    case 9:
        ui->postPVar2CBox->addItem("Porosity","MatPorosity");
        ui->postPVar2CBox->addItem("E","MatE");
        ui->postPVar2CBox->addItem("Poisson","MatPoisson");
    break;

    default:
    break;
    } //switch

    //qDebug() << "Saiu populate_postPVar2CBox";
}

void MainWindow::changedTabs_updateCheckboxes (bool sts) {
    // Get list of tabs (dockwidgets) stacked together
    QList<QTabBar *> tabList = this->findChildren<QTabBar *>();
    QTabBar *tabBar = tabList.at(1); // Checar e talvez mudar quando adicionar novas tablists
    int idx_local = tabBar->currentIndex();

#ifdef DEBUG
    qDebug() << "#######################################################" << tabBar->tabText(idx_local);
#endif

    changedTabs_updateCheckboxes(idx_local);
}

void MainWindow::changedTabs_updateCheckboxes (int idx) {
#ifdef DEBUG
    qDebug() << "Entrou changedTabs_updateCheckboxes idx=" << idx;
#endif

    // PWB tab
    switch (idx) {
        case 0: //PWB
            ellipsCheck_toggled(false);
            polyChainCheck_toggled(false);
            contourCheck_toggled(false);
        break;
        case 1: //Breakout
            ellipsCheck_toggled(identifyellipsDock->ui->ellipsCheck->isChecked());
            polyChainCheck_toggled(identifyellipsDock->ui->polyChainCheck->isChecked());
            contourCheck_toggled(identifyellipsDock->ui->contourCheck->isChecked());
        break;
        case 2: //P Ref
            ellipsCheck_toggled(false);
            polyChainCheck_toggled(prefineDock->ui->polyChainCheck->isChecked());
            contourCheck_toggled(prefineDock->ui->contourCheck->isChecked());
        break;
        case 3: //H Ref
            ellipsCheck_toggled(false);
            polyChainCheck_toggled(hrefineDock->ui->polyChainCheck->isChecked());
            contourCheck_toggled(hrefineDock->ui->contourCheck->isChecked());
        break;
        case 4: //Acid job
            ellipsCheck_toggled(false);
            polyChainCheck_toggled(false);
            contourCheck_toggled(false);
        break;
        case 5: //Casing
            ellipsCheck_toggled(false);
            polyChainCheck_toggled(false);
            contourCheck_toggled(false);
        break;
        case 6: //Solution
            ellipsCheck_toggled(false);
            polyChainCheck_toggled(false);
            contourCheck_toggled(false);
        break;
        default: // do nothing
        break;
    }
#ifdef DEBUG
    qDebug() << "Saiu changedTabs_updateCheckboxes";
#endif

}

void MainWindow::ShowHistoryContextMenu(const QPoint &pos)
{
#ifdef DEBUG
    qDebug() << "Entrou contexto";
#endif
    QListWidgetItem *Item = historyDock->ui->historyList->itemAt(pos);
    if (!Item) return;

    // se eh um sub item, nao mostrar menu
    if (Item->text().toStdString().c_str()[0] == ' ') return;

    QPoint globalPos = historyDock->ui->historyList->mapToGlobal(pos); // for most widgets

    //Create context menu and its actions
    QMenu myMenu;
    QAction *aReRunFromHere = myMenu.addAction(tr("Rollback to here."));

    // get selected action
    QAction* selectedItem = myMenu.exec(globalPos);

    if (aReRunFromHere == selectedItem) {
        int ret_idx = Item->data(5).toInt();
#ifdef DEBUG
            qDebug() << "Escolheu rollback to id = " << ret_idx << " list size = " << fwellb.GetConfigListSize();
#endif
        while (ret_idx < fwellb.GetConfigListSize()-1) {
            fwellb.PopConfiguration();
        }
    }

    fLoadedConfig = fwellb.GetCurrentConfig();

#ifdef DEBUG
    qDebug() << "AKIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII";
#endif

    //composeScene();

    openBin();

    populate_GUIFields();

    reloadHistory();

#ifdef DEBUG
    qDebug() << "Loaded id = (list size - 1) = " << fwellb.GetConfigListSize() -1;
#endif

#ifdef DEBUG
    qDebug() << "Saiu contexto";
#endif
}