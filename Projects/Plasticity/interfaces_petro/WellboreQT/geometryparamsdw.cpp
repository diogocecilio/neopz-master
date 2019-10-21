#include "geometryparamsdw.h"
#include "ui_geometryparamsdw.h"

#ifdef DEBUG
#include <QDebug>
#endif

geometryparamsDW::geometryparamsDW(QWidget *parent) :
    QDockWidget(parent),
    ui(new Ui::geometryparamsDW)
{
    ui->setupUi(this);

    fWB_Dir = EVerticalWell;

    connect (ui->rb_direction_SH, SIGNAL(toggled(bool)), this, SLOT(wellboreDirectionChanged()));
    connect (ui->rb_direction_Sh, SIGNAL(toggled(bool)), this, SLOT(wellboreDirectionChanged()));
    connect (ui->rb_direction_SV, SIGNAL(toggled(bool)), this, SLOT(wellboreDirectionChanged()));

    //update image
    setDirection(0);
    wellboreDirectionChanged();
}

geometryparamsDW::~geometryparamsDW()
{
    delete ui;
}

void geometryparamsDW::wellboreDirectionChanged() {
    if (ui->rb_direction_SV->isChecked()) {
        // TODO
        fWB_Dir = EVerticalWell;
        ui->wellbore_direction_img->setText("<img src=':imgs/wellsigmav.png' height='200' />");
    }

    if (ui->rb_direction_SH->isChecked()) {
        // TODO
        fWB_Dir = EHorizontalWellalongH;
        ui->wellbore_direction_img->setText("<img src=':imgs/wellsigmaH1.png' height='200' />");
    }

    if (ui->rb_direction_Sh->isChecked()) {
        // TODO
        fWB_Dir = EHorizontalWellalongh;
        ui->wellbore_direction_img->setText("<img src=':imgs/wellsigmah.png' height='200' />");
    }

#ifdef DEBUG
    qDebug() << "New Wellbore direction value: " << fWB_Dir;
#endif
}

void geometryparamsDW::setDirection(int conf) {

    if (conf == EVerticalWell) {
        ui->rb_direction_SV->isChecked();
        fWB_Dir = EVerticalWell;
    }

    if (conf == EHorizontalWellalongH) {
        ui->rb_direction_SH->isChecked();
        fWB_Dir = EHorizontalWellalongH;
    }

    if (conf == EHorizontalWellalongh) {
        ui->rb_direction_Sh->isChecked();
        fWB_Dir = EHorizontalWellalongh;
    }

#ifdef DEBUG
    qDebug() << "New Wellbore direction value: " << fWB_Dir;
#endif
}

int geometryparamsDW::getDirection() {
    return fWB_Dir;
}