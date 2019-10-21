#ifndef MESHDW_H
#define MESHDW_H

#include <QDockWidget>

namespace Ui {
class meshDW;
}

class meshDW : public QDockWidget
{
    Q_OBJECT
    
public:
    explicit meshDW(QWidget *parent = 0);
    ~meshDW();
    
//private:
    Ui::meshDW *ui;
};

#endif // MESHDW_H
