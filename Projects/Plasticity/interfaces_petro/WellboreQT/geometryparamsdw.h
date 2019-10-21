#ifndef GEOMETRYPARAMSDW_H
#define GEOMETRYPARAMSDW_H

#include <QDockWidget>

namespace Ui {
class geometryparamsDW;
}

class geometryparamsDW : public QDockWidget
{
    Q_OBJECT

    /// Indication of the well configuration
    enum EWellConfiguration {ENoConfig, EVerticalWell, EHorizontalWellalongh, EHorizontalWellalongH};

    EWellConfiguration fWB_Dir;

private slots:
    void wellboreDirectionChanged();
    
public:
    explicit geometryparamsDW(QWidget *parent = 0);
    ~geometryparamsDW();
    void setDirection(int idx);
    int getDirection();
    
//private:
    Ui::geometryparamsDW *ui;
};

#endif // GEOMETRYPARAMSDW_H
