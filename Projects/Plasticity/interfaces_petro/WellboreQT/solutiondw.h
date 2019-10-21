#ifndef SOLUTIONDW_H
#define SOLUTIONDW_H

#include <QDockWidget>

namespace Ui {
class solutionDW;
}

class solutionDW : public QDockWidget
{
    Q_OBJECT

public:
    explicit solutionDW(QWidget *parent = 0);
    ~solutionDW();

//private:
    Ui::solutionDW *ui;
};

#endif // SOLUTIONDW_H
