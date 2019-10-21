#include "fluiddw.h"
#include "ui_fluiddw.h"

#include "greekletters.h"

fluidDW::fluidDW(QWidget *parent) :
    QDockWidget(parent),
    ui(new Ui::fluidDW)
{
    ui->setupUi(this);

    //ui->labelViscosity->setText(QString("Viscosity ").append(mu).append("="));
    ui->labelViscosity->setText(ui->labelViscosity->text().replace("mu", mu));
}

fluidDW::~fluidDW()
{
    delete ui;
}
