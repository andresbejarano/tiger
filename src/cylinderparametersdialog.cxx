#include "cylinderparametersdialog.h"
#include "ui_cylinderparametersdialog.h"

CylinderParametersDialog::CylinderParametersDialog(QWidget *parent) :
    QDialog(parent),
    ui(new Ui::CylinderParametersDialog)
{
    ui->setupUi(this);
}

CylinderParametersDialog::~CylinderParametersDialog()
{
    delete ui;
}

QString CylinderParametersDialog::GetAxis() const
{
    return ui->axisComboBox->currentText();
}

double CylinderParametersDialog::GetLength() const
{
    return ui->lengthDoubleSpinBox->value();
}

size_t CylinderParametersDialog::GetLengthSegments() const
{
    return ui->lengthSegmentsSpinBox->value();
}

size_t CylinderParametersDialog::GetRadialSegments() const
{
    return ui->radialSegmentsSpinBox->value();
}

double CylinderParametersDialog::GetRadius() const
{
    return ui->radiusDoubleSpinBox->value();
}
