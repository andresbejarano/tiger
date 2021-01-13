#include "coneparametersdialog.h"
#include "ui_coneparametersdialog.h"

ConeParametersDialog::ConeParametersDialog(QWidget *parent) :
    QDialog(parent),
    ui(new Ui::ConeParametersDialog)
{
    ui->setupUi(this);
}

ConeParametersDialog::~ConeParametersDialog()
{
    delete ui;
}

QString ConeParametersDialog::GetAxis() const
{
    return ui->axisComboBox->currentText();
}

double ConeParametersDialog::GetLength() const
{
    return ui->lengthDoubleSpinBox->value();
}

size_t ConeParametersDialog::GetLengthSegments() const
{
    return ui->lengthSegmentsSpinBox->value();
}

size_t ConeParametersDialog::GetRadialSegments() const
{
    return ui->radialSegmentsSpinBox->value();
}

double ConeParametersDialog::GetRadius() const
{
    return ui->radiusDoubleSpinBox->value();
}
