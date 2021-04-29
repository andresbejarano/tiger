#include "squaregridparametersdialog.h"
#include "ui_squaregridparametersdialog.h"

SquareGridParametersDialog::SquareGridParametersDialog(QWidget *parent) :
    QDialog(parent),
    ui(new Ui::SquareGridParametersDialog)
{
    ui->setupUi(this);
}

SquareGridParametersDialog::~SquareGridParametersDialog()
{
    delete ui;
}

double SquareGridParametersDialog::GetWidth() const
{
    return ui->widthDoubleSpinBox->value();
}

double SquareGridParametersDialog::GetHeight() const
{
    return ui->heightDoubleSpinBox->value();
}

size_t SquareGridParametersDialog::GetWidthSegments() const
{
    return ui->widthSegmentsSpinBox->value();
}

size_t SquareGridParametersDialog::GetHeightSegments() const
{
    return ui->heightSegmentsSpinBox->value();
}

QString SquareGridParametersDialog::GetPlane() const
{
    return ui->planeComboBox->currentText();
}
