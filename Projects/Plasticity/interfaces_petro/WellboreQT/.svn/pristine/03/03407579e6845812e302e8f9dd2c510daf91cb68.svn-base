#include "casingdw.h"
#include "ui_casingdw.h"

#include "greekletters.h"

#ifdef DEBUG
#include <QDebug>
#endif

casingDW::casingDW(QWidget *parent) :
    QDockWidget(parent),
    ui(new Ui::casingDW)
{
    ui->setupUi(this);

    // Linear elastic material
    ui->nu_label->setText(ui->nu_label->text().replace("nu", nu));

    ui->vonMises_GBox->hide();
    fMat = 0; //linear = default

    connect (ui->material_CBox, SIGNAL(currentIndexChanged(int)), this, SLOT(ChangeMaterial(int)));
}

void casingDW::ChangeMaterial(int idx) {
#ifdef DEBUG
    qDebug() << "ChangeMaterialType idx = " << idx;
#endif
    if (idx == 0) {
        ui->vonMises_GBox->hide();
        ui->linearElastic_GBox->show();
#ifdef DEBUG
        qDebug() << "linear elastic";
#endif
    }
    if (idx == 1) {
        ui->linearElastic_GBox->hide();
        ui->vonMises_GBox->show();
#ifdef DEBUG
        qDebug() << "von mises";
#endif
    }
    fMat = idx;
}

int casingDW::SelectedMat() {
    return fMat;
}

void casingDW::selectMat(int idx) {
    ui->material_CBox->setCurrentIndex(idx);
    ChangeMaterial(idx);
}

casingDW::~casingDW()
{
    delete ui;
}
