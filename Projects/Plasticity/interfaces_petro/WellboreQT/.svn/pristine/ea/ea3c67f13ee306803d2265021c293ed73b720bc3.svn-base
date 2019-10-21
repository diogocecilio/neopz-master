#include "initialconfigdialog.h"
#include "ui_initialconfigdialog.h"

initialConfigDialog::initialConfigDialog(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::initialConfigDialog)
{
    this ->setWindowModality(Qt::ApplicationModal);

    ui->setupUi(this);

    connect (this->ui->cancelBtn, SIGNAL(clicked()), this, SLOT(close()) );
}

initialConfigDialog::~initialConfigDialog()
{
    delete ui;
}
