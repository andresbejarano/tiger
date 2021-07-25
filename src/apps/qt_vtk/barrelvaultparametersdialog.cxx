#include "barrelvaultparametersdialog.h"
#include "ui_barrelvaultparametersdialog.h"

BarrelVaultParametersDialog::BarrelVaultParametersDialog(QWidget *parent) :
    QDialog(parent),
    ui(new Ui::BarrelVaultParametersDialog)
{
    ui->setupUi(this);
}

BarrelVaultParametersDialog::~BarrelVaultParametersDialog()
{
    delete ui;
}

double BarrelVaultParametersDialog::GetLength() const
{
    return ui->lengthDoubleSpinBox->value();
}

double BarrelVaultParametersDialog::GetRadius() const
{
    return ui->radiusDoubleSpinBox->value();
}

size_t BarrelVaultParametersDialog::GetLengthSegments() const
{
    return ui->lengthSegmentsSpinBox->value();
}

size_t BarrelVaultParametersDialog::GetRadialSegments() const
{
    return ui->radialSegmentsSpinBox->value();
}

QString BarrelVaultParametersDialog::GetAxis() const
{
    return ui->axisComboBox->currentText();
}
