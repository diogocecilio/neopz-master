#ifndef PWBDW_H
#define PWBDW_H

#include <QDockWidget>

namespace Ui {
class pwbDW;
}

class pwbDW : public QDockWidget
{
    Q_OBJECT

private slots:
    void fluidModelChanged();

signals:
    void fluidModelChanged(int FluidModel);
    
public:
    explicit pwbDW(QWidget *parent = 0);
    ~pwbDW();
    
//private:
    Ui::pwbDW *ui;

private:
    int fFluidModel;
};

#endif // PWBDW_H
