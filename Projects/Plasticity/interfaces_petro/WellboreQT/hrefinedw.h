#ifndef HREFINEDW_H
#define HREFINEDW_H

#include <QDockWidget>

namespace Ui {
class hrefineDW;
}

class hrefineDW : public QDockWidget
{
    Q_OBJECT
    
public:
    explicit hrefineDW(QWidget *parent = 0);
    ~hrefineDW();
    
//private:
    Ui::hrefineDW *ui;
};

#endif // HREFINEDW_H
