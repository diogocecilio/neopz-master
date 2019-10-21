#ifndef POROMECHDW_H
#define POROMECHDW_H

#include <QDockWidget>

namespace Ui {
class poromechDW;
}

class poromechDW : public QDockWidget
{
    Q_OBJECT

private slots:
    void recalculateA();

public:
    explicit poromechDW(QWidget *parent = 0);
    ~poromechDW();

    double getA () {
        return fA;
    }
    
//private:
    Ui::poromechDW *ui;

private:
    double fA;
};

#endif // POROMECHDW_H
