#include "truncatedconeparametersdialog.h"
#include "ui_truncatedconeparametersdialog.h"

TruncatedConeParametersDialog::TruncatedConeParametersDialog(QWidget *parent) :
    QDialog(parent),
    ui(new Ui::TruncatedConeParametersDialog)
{
    ui->setupUi(this);
}

TruncatedConeParametersDialog::~TruncatedConeParametersDialog()
{
    delete ui;
}

QString TruncatedConeParametersDialog::GetAxis() const
{
    return ui->axisComboBox->currentText();
}

double TruncatedConeParametersDialog::GetBottomRadius() const
{
    return ui->bottomRadiusDoubleSpinBox->value();
}

double TruncatedConeParametersDialog::GetLength() const
{
    return ui->lengthDoubleSpinBox->value();
}

size_t TruncatedConeParametersDialog::GetLengthSegments() const
{
    return ui->lengthSegmentsSpinBox->value();
}

size_t TruncatedConeParametersDialog::GetRadialSegments() const
{
    return ui->radialSegmentsSpinBox->value();
}

double TruncatedConeParametersDialog::GetTopRadius() const
{
    return ui->topRadiusDoubleSpinBox->value();
}
