#include "paraboloidparametersdialog.h"
#include "ui_paraboloidparametersdialog.h"

ParaboloidParametersDialog::ParaboloidParametersDialog(QWidget *parent) :
    QDialog(parent),
    ui(new Ui::ParaboloidParametersDialog)
{
    ui->setupUi(this);
}

ParaboloidParametersDialog::~ParaboloidParametersDialog()
{
    delete ui;
}

double ParaboloidParametersDialog::GetWidth() const
{
    return ui->widthDoubleSpinBox->value();
}

double ParaboloidParametersDialog::GetHeight() const
{
    return ui->heightDoubleSpinBox->value();
}

size_t ParaboloidParametersDialog::GetWidthSegments() const
{
    return ui->widthSegmentsSpinBox->value();
}

size_t ParaboloidParametersDialog::GetHeightSegments() const
{
    return ui->heightSegmentsSpinBox->value();
}

QString ParaboloidParametersDialog::GetPlane() const
{
    return ui->planeComboBox->currentText();
}

double ParaboloidParametersDialog::GetA() const
{
    return ui->aDoubleSpinBox->value();
}

double ParaboloidParametersDialog::GetB() const
{
    return ui->bDoubleSpinBox->value();
}

QString ParaboloidParametersDialog::GetType() const
{
    return ui->typeComboBox->currentText();
}
