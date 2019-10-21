#include "prefinedw.h"
#include "ui_prefinedw.h"

prefineDW::prefineDW(QWidget *parent) :
    QDockWidget(parent),
    ui(new Ui::prefineDW)
{
    ui->setupUi(this);
}

prefineDW::~prefineDW()
{
    delete ui;
}
