#ifndef PARAMETERSDW_H
#define PARAMETERSDW_H

#include <QDockWidget>

namespace Ui {
class parametersDW;
}

class parametersDW : public QDockWidget
{
    Q_OBJECT

    /// holds the idx of selected material
    ///  0=Sandler, 1=MohrCoulomb, 2=Elastic
    int fMat;

private slots:
    void ChangeMaterial(int idx);

public:
    explicit parametersDW(QWidget *parent = 0);
    ~parametersDW();
    int SelectedMat();
    void selectMat(int idx);

//private:
    Ui::parametersDW *ui;
};

#endif // PARAMETERSDW_H
