#include "evaluatedw.h"
#include "ui_evaluatedw.h"

evaluateDW::evaluateDW(QWidget *parent) :
    QDockWidget(parent),
    ui(new Ui::evaluateDW)
{
    ui->setupUi(this);
}

evaluateDW::~evaluateDW()
{
    delete ui;
}
