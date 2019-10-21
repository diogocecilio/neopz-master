#ifndef INITIALCONFIGDIALOG_H
#define INITIALCONFIGDIALOG_H

#include <QMainWindow>

namespace Ui {
class initialConfigDialog;
}

class initialConfigDialog : public QMainWindow
{
    Q_OBJECT
    
public:
    explicit initialConfigDialog(QWidget *parent = 0);
    ~initialConfigDialog();

//private slots:
    //materialTypeChanged();
    
//private:
    Ui::initialConfigDialog *ui;
};

#endif // INITIALCONFIGDIALOG_H
