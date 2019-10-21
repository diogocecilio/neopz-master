#include "pwbdw.h"
#include "ui_pwbdw.h"

#include "WellBoreAnalysis.h"

pwbDW::pwbDW(QWidget *parent) :
    QDockWidget(parent),
    ui(new Ui::pwbDW)
{
    ui->setupUi(this);

    fFluidModel = ENonPenetrating;

    connect (ui->FluidModelRB_P, SIGNAL(toggled(bool)), this, SLOT(fluidModelChanged()));
    connect (ui->FluidModelRB_NP, SIGNAL(toggled(bool)), this, SLOT(fluidModelChanged()));
    connect (ui->FluidModelRB_Const, SIGNAL(toggled(bool)), this, SLOT(fluidModelChanged()));
}

pwbDW::~pwbDW()
{
    delete ui;
}

void pwbDW::fluidModelChanged() {
    if (ui->FluidModelRB_P->isChecked()) {
        if (fFluidModel != EPenetrating) {
            fFluidModel = EPenetrating;
            emit fluidModelChanged(EPenetrating);
            return;
        }
    }
    if (ui->FluidModelRB_NP->isChecked()) {
        if (fFluidModel != ENonPenetrating) {
            fFluidModel = ENonPenetrating;
            emit fluidModelChanged(ENonPenetrating);
            return;
        }
    }
    if (ui->FluidModelRB_Const->isChecked()) {
        //TODO
        return;
    }
}
