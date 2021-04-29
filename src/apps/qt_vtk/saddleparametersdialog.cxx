#include "saddleparametersdialog.h"
#include "ui_saddleparametersdialog.h"

SaddleParametersDialog::SaddleParametersDialog(QWidget *parent) :
    QDialog(parent),
    ui(new Ui::SaddleParametersDialog)
{
    ui->setupUi(this);
}

SaddleParametersDialog::~SaddleParametersDialog()
{
    delete ui;
}

double SaddleParametersDialog::GetWidth() const
{
    return ui->widthDoubleSpinBox->value();
}

double SaddleParametersDialog::GetHeight() const
{
    return ui->heightDoubleSpinBox->value();
}

size_t SaddleParametersDialog::GetWidthSegments() const
{
    return ui->widthSegmentsSpinBox->value();
}

size_t SaddleParametersDialog::GetHeightSegments() const
{
    return ui->heightSegmentsSpinBox->value();
}

QString SaddleParametersDialog::GetPlane() const
{
    return ui->planeComboBox->currentText();
}
