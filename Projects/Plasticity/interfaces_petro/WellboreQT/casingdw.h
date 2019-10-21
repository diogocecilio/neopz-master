#ifndef CASINGDW_H
#define CASINGDW_H

#include <QDockWidget>

namespace Ui {
class casingDW;
}

class casingDW : public QDockWidget
{
    Q_OBJECT

    /// holds the idx of selected material
    ///  0=Elastic , 1=VonMises
    int fMat;

private slots:
    void ChangeMaterial(int idx);

public:
    explicit casingDW(QWidget *parent = 0);
    ~casingDW();
    int SelectedMat();
    void selectMat(int idx);

//private:
    Ui::casingDW *ui;
};

#endif // CASINGDW_H
