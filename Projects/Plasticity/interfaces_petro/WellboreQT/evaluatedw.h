#ifndef EVALUATEDW_H
#define EVALUATEDW_H

#include <QDockWidget>

namespace Ui {
class evaluateDW;
}

class evaluateDW : public QDockWidget
{
    Q_OBJECT
    
public:
    explicit evaluateDW(QWidget *parent = 0);
    ~evaluateDW();
    
//private:
    Ui::evaluateDW *ui;
};

#endif // EVALUATEDW_H
