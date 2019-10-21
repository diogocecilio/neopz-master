#include "tpbrunitinput.h"

#include <QGridLayout>
#include <QLineEdit>
#include <QLabel>
#include <QComboBox>
#include <QHBoxLayout>
#include <QVariant>
//#include <QDebug>
extern "C" {
    #include "externallibs/units-2.11/units.h"
}

// Units fix: create a temporary file with contents of :/definition.units (resource file)
QString TPBrUnitInput::fTmpDefsUnitsFile;

TPBrUnitInput::TPBrUnitInput(QWidget *parent) :
    QWidget(parent)
{
    fSIValue = 0.;
    fTVD = 0.;

    fPrecision = 3;

    fDisplayValue = new QLineEdit;
    fDisplayValue->setText(QVariant(fSIValue).toString());

    fDisplayUnit = new QComboBox;

    setUnitType(EDimensionless);

    QHBoxLayout *mainLayout = new QHBoxLayout;
    mainLayout->addWidget(fDisplayValue);
    mainLayout->addWidget(fDisplayUnit);
    mainLayout->addStretch();

    setLayout(mainLayout);

    connect(fDisplayValue, SIGNAL(editingFinished()), this, SLOT(set()));
    connect(fDisplayUnit, SIGNAL(currentIndexChanged(int)), this, SLOT(convert(int)));
}

TPBrUnitInput::~TPBrUnitInput () {
    delete fDisplayValue;
    delete fDisplayUnit;
    delete this->layout();
}

void TPBrUnitInput::reload(double newTVD) {
    fTVD = newTVD;
    //Erick asked to not reflect the new TVD in displayed value
    //convert(fDisplayUnit->currentIndex());
}

void TPBrUnitInput::changeUnit(int idx) {
    fDisplayUnit->setCurrentIndex(idx);
}

// Display converted internal value according to combobox option
void TPBrUnitInput::convert(int idx) {
    if (fUnitType != EDimensionless) {
        QString fromStr(fDisplayUnit->itemData(idx).toString());
        QString toStr(fDisplayUnit->itemText(idx));

        double multiplier = 0.;
        double divider = 0.;

        // All exceptions
        if (toStr == "ppg") {
            // Pre-conversion: Convert internal value to psi
            int rett = convertGNUUnits(fTmpDefsUnitsFile.toStdString().c_str(), fromStr.toStdString().c_str(), "psi", &multiplier, &divider);
            //P[psi] = 0,170433233 * Grad[ppg] * Prof [m]
            //=> Grad[ppg] = P[psi] / (0,170433233 * Prof [m])
            multiplier = multiplier / (0.170433233 * fTVD);

            fDisplayValue->setText(QString::number( fSIValue * multiplier, 'f', fPrecision ));
            return;
        }

        int rett = convertGNUUnits(fTmpDefsUnitsFile.toStdString().c_str(), fromStr.toStdString().c_str(), toStr.toStdString().c_str(), &multiplier, &divider);

        if (!rett) {
            fDisplayValue->setText(QString::number( fSIValue * multiplier, 'f', fPrecision ));
        }
    }
    else
    {
        fDisplayValue->setText(QString::number( fSIValue, 'f', fPrecision ));
    }
}

// Change internal value, converting it according to combobox option
void TPBrUnitInput::set() {
    if (fUnitType != EDimensionless) {
        QString fromStr(fDisplayUnit->itemText(fDisplayUnit->currentIndex()));
        QString toStr(fDisplayUnit->itemData(fDisplayUnit->currentIndex()).toString());

        double multiplier = 0.;
        double divider = 0.;

        // All exceptions
        if (fromStr == "ppg") {
            //P[psi] = 0,170433233 * Grad[ppg] * Prof [m]
            double tmpPsi = 0.170433233 * fDisplayValue->text().toDouble() * fTVD;
            // tmp conversion: Convert temporary psi value to internal value
            int rett = convertGNUUnits(fTmpDefsUnitsFile.toStdString().c_str(), "psi", "MPa", &multiplier, &divider);
            fSIValue = QString::number( tmpPsi * multiplier, 'f', fPrecision ).toDouble();

            emit valueChanged(fSIValue);
            return;
        }

        int rett = convertGNUUnits(fTmpDefsUnitsFile.toStdString().c_str(), fromStr.toStdString().c_str(), toStr.toStdString().c_str(), &multiplier, &divider);

        if (!rett) {
            fSIValue = QString::number( fDisplayValue->text().toDouble() * multiplier, 'f', fPrecision ).toDouble();
        }
    }
    else
    {
        fSIValue = QString::number( fDisplayValue->text().toDouble(), 'f', fPrecision ).toDouble();
    }

    emit valueChanged(fSIValue);
}

QString TPBrUnitInput::text() const {
    return QString::number( fSIValue, 'f', fPrecision );
}

void TPBrUnitInput::setText( const QString &text ) {
    fSIValue = QString::number( text.toDouble(), 'f', fPrecision ).toDouble();
    fDisplayValue->setText(QString::number( fSIValue, 'f', fPrecision ));

    emit valueChanged(fSIValue);
}

QString TPBrUnitInput::currentUnitText() const {
    return fDisplayUnit->itemData(fDisplayUnit->currentIndex()).toString();
}

TPBrUnitInput::EUnitType TPBrUnitInput::unitType() const {
    return fUnitType;
}

void TPBrUnitInput::setUnitType(const EUnitType & unittp) {

    //if (unittp == fUnitType) return;

    fUnitType = unittp;

    fDisplayUnit->clear();
    /*
    fDisplayUnit->addItem("FROM UNIT", QVariant("CONVERT TO 'SI' UNIT"));
    */
    if (fUnitType == EPressure) {
        fDisplayUnit->addItem("MPa", QVariant("MPa"));
        fDisplayUnit->addItem("Pa", QVariant("MPa"));
        fDisplayUnit->addItem("psi", QVariant("MPa"));
        fDisplayUnit->addItem("ppg", QVariant("MPa"));
    }
    if (fUnitType == EDimension) {
        fDisplayUnit->addItem("m", QVariant("m"));
        fDisplayUnit->addItem("ft", QVariant("m"));
        fDisplayUnit->addItem("in", QVariant("m"));
    }
    if (fUnitType == EAngle) {
        fDisplayUnit->addItem("radian", QVariant("radian"));
        fDisplayUnit->addItem("degree", QVariant("radian"));
    }
    if (fUnitType == EDimensionless) {
        fDisplayUnit->hide();
        return;
    }
    fDisplayUnit->show();
}

bool TPBrUnitInput::isEnabled() const {
    return fDisplayValue->isEnabled();
}

void TPBrUnitInput::setEnabled ( bool status ) {
    fDisplayValue->setEnabled(status);
    fDisplayUnit->setEnabled(status);
}
