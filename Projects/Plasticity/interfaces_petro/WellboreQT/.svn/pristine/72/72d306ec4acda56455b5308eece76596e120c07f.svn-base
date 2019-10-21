#include "acidjobdw.h"
#include "ui_acidjobdw.h"

#include "greekletters.h"

acidjobDW::acidjobDW(QWidget *parent) :
    QDockWidget(parent),
    ui(new Ui::acidjobDW)
{
    ui->setupUi(this);

    // Weakening tab
    ui->E_formula_label->setText(ui->E_formula_label->text().replace("Phi", Phi).replace("Delta", Delta));

    // Rock tab
    ui->Phi0_label->setText(ui->Phi0_label->text().replace("Phi", Phi));
    ui->Phimax_label->setText(ui->Phimax_label->text().replace("Phi", Phi));

    // Reaction tab
    ui->betaf_label->setText(ui->betaf_label->text().replace("beta", beta));
    ui->betam_label->setText(ui->betam_label->text().replace("beta", beta));
    ui->rhom_label->setText(ui->rhom_label->text().replace("rho", rho));

}

acidjobDW::~acidjobDW()
{
    delete ui;
}
