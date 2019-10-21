#ifndef TPBRUNITINPUT_H
#define TPBRUNITINPUT_H

#include <QWidget>
#include <QTemporaryFile>

class QLineEdit;
class QLabel;
class QComboBox;
class QString;

class TPBrUnitInput : public QWidget
{
    Q_OBJECT

    Q_ENUMS (EUnitType)

    Q_PROPERTY( QString text READ text WRITE setText )
    Q_PROPERTY( EUnitType unitType READ unitType WRITE setUnitType )
    Q_PROPERTY(bool enabled READ isEnabled WRITE setEnabled)

public:

    enum EUnitType {EPressure, EDimension, EAngle, EDimensionless};

    explicit TPBrUnitInput(QWidget *parent = 0);
    virtual ~TPBrUnitInput();

    QString text() const;
    void setText( const QString & text);

    QString currentUnitText() const;

    EUnitType unitType() const;
    void setUnitType(const EUnitType & unittp);

    bool isEnabled() const;
    void setEnabled ( bool status );

    // Units fix: create a temporary file with contents of :/definition.units (resource file)
    static QString fTmpDefsUnitsFile;
    
signals:
    void valueChanged(double newValue);
    
public Q_SLOTS:
    void reload(double);
    void changeUnit(int idx);

private:
    QLineEdit *fDisplayValue;
    QComboBox *fDisplayUnit;

    double fSIValue;
    EUnitType fUnitType;

    double fTVD;

    int fPrecision;

private slots:
    void convert(int idx);
    void set();
    
};

#endif // TPBRUNITINPUT_H
