#include "poromechdw.h"
#include "ui_poromechdw.h"

#include "greekletters.h"

#ifdef DEBUG
#include <QDebug>
#endif

poromechDW::poromechDW(QWidget *parent) :
    QDockWidget(parent),
    ui(new Ui::poromechDW)
{
    ui->setupUi(this);

    fA = 0.;

    ui->labelAlphaFormula->setText(QString("%1 (1-2%2) / (1-%3)").arg(alpha, nu, nu));
    //ui->labelBiotWillis->setText(QString("off-Boundary Biot-Willis %1:").arg(alpha));
    ui->labelBiotWillis->setText(ui->labelBiotWillis->text().replace("alpha", alpha));

    connect(ui->alpha, SIGNAL(editingFinished()), this, SLOT(recalculateA()));
    connect(ui->nu, SIGNAL(editingFinished()), this, SLOT(recalculateA()));
    connect(ui->A, SIGNAL(editingFinished()), this, SLOT(recalculateA()));
    connect(ui->aFormulaRadio, SIGNAL(toggled(bool)), this, SLOT(recalculateA()));
    connect(ui->aUserRadio, SIGNAL(toggled(bool)), this, SLOT(recalculateA()));
    connect(ui->aBCRadio, SIGNAL(toggled(bool)), this, SLOT(recalculateA()));

    connect(ui->aUserRadio, SIGNAL(toggled(bool)), ui->alpha, SLOT(setDisabled(bool)));
    connect(ui->aUserRadio, SIGNAL(toggled(bool)), ui->nu, SLOT(setDisabled(bool)));

    connect(ui->aFormulaRadio, SIGNAL(toggled(bool)), ui->A, SLOT(setDisabled(bool)));

    ui->A->setDisabled(true);
    recalculateA();
}

poromechDW::~poromechDW()
{
    delete ui;
}

void poromechDW::recalculateA() {

    if (ui->aFormulaRadio->isChecked()) {
        double A = 0;

        double alpha = ui->alpha->text().toDouble();
        double nu = ui->nu->text().toDouble();

        A = alpha * ((1-2*nu) / (1-nu));

        fA = A;
        ui->A->setText(QVariant(A).toString());
    }

    if (ui->aUserRadio->isChecked()) {
        fA = ui->A->text().toDouble();
    }

    if (ui->aBCRadio->isChecked()) {
        // TODO
    }

#ifdef DEBUG
    qDebug() << "New A value: " << fA;
#endif

}
