#include "hrefinedw.h"
#include "ui_hrefinedw.h"

hrefineDW::hrefineDW(QWidget *parent) :
    QDockWidget(parent),
    ui(new Ui::hrefineDW)
{
    ui->setupUi(this);
}

hrefineDW::~hrefineDW()
{
    delete ui;
}
