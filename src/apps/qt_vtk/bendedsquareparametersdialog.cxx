#include "bendedsquareparametersdialog.h"
#include "ui_bendedsquareparametersdialog.h"

BendedSquareParametersDialog::BendedSquareParametersDialog(QWidget *parent) :
    QDialog(parent),
    ui(new Ui::BendedSquareParametersDialog)
{
    ui->setupUi(this);
}

BendedSquareParametersDialog::~BendedSquareParametersDialog()
{
    delete ui;
}

QString BendedSquareParametersDialog::GetBendDirection() const
{
    return ui->bendDirectionComboBox->currentText();
}

double BendedSquareParametersDialog::GetLength() const
{
    return ui->lengthDoubleSpinBox->value();
}

size_t BendedSquareParametersDialog::GetLengthSegments() const
{
    return ui->lengthSegmentsSpinBox->value();
}

size_t BendedSquareParametersDialog::GetRadialSegments() const
{
    return ui->radialSegmentsSpinBox->value();
}

double BendedSquareParametersDialog::GetRadius() const
{
    return ui->radiusDoubleSpinBox->value();
}
