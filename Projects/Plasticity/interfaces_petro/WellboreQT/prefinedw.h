#ifndef PREFINEDW_H
#define PREFINEDW_H

#include <QDockWidget>

namespace Ui {
class prefineDW;
}

class prefineDW : public QDockWidget
{
    Q_OBJECT
    
public:
    explicit prefineDW(QWidget *parent = 0);
    ~prefineDW();
    
//private:
    Ui::prefineDW *ui;
};

#endif // PREFINEDW_H
