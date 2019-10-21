#include "meshdw.h"
#include "ui_meshdw.h"

#include "greekletters.h"

meshDW::meshDW(QWidget *parent) :
    QDockWidget(parent),
    ui(new Ui::meshDW)
{
    ui->setupUi(this);

    //ui->labelFirstElem->setText(QString("").append(Delta).append("r /").append(Delta).append("c 1st element row:"));
    ui->labelFirstElem->setText(ui->labelFirstElem->text().replace("delta", Delta));
}

meshDW::~meshDW()
{
    delete ui;
}
