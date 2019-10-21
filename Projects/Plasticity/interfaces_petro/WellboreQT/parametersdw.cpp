#include "parametersdw.h"
#include "ui_parametersdw.h"

#include "greekletters.h"

#ifdef DEBUG
#include <QDebug>
#endif

parametersDW::parametersDW(QWidget *parent) :
    QDockWidget(parent),
    ui(new Ui::parametersDW)
{
    ui->setupUi(this);

    ui->mohrCGrp->hide();
    ui->elasticGrp->hide();
    fMat = 0; //sandler = default

    ui->materialCBox->setEnabled(true);

    //ui->labelBiotWillis->setText(QString("Biot-Willis %1:").arg(alpha));
    ui->labelBiotWillis->setText(ui->labelBiotWillis->text().replace("a", alpha));

    connect (ui->materialCBox, SIGNAL(currentIndexChanged(int)), this, SLOT(ChangeMaterial(int)));
}

void parametersDW::ChangeMaterial(int idx) {
#ifdef DEBUG
    qDebug() << "ChangeMaterialType idx = " << idx;
#endif
    if (idx == 0) {
        ui->mohrCGrp->hide();
        ui->elasticGrp->hide();
        ui->sandlerGrp->show();
#ifdef DEBUG
        qDebug() << "Sandler";
#endif
    }
    if (idx == 1) {
        ui->sandlerGrp->hide();
        ui->elasticGrp->hide();
        ui->mohrCGrp->show();
#ifdef DEBUG
        qDebug() << "mohr";
#endif
    }
    if (idx == 2) {
        ui->sandlerGrp->hide();
        ui->mohrCGrp->hide();
        ui->elasticGrp->show();
#ifdef DEBUG
        qDebug() << "elastico";
#endif
    }
    fMat = idx;
}

int parametersDW::SelectedMat() {
    return fMat;
}

void parametersDW::selectMat(int idx) {
    ui->materialCBox->setCurrentIndex(idx);
    ChangeMaterial(idx);
}

parametersDW::~parametersDW()
{
    delete ui;
}
