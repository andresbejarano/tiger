#include "equilateraltriangleparametersdialog.h"
#include "ui_equilateraltriangleparametersdialog.h"

EquilateralTriangleParametersDialog::EquilateralTriangleParametersDialog(QWidget *parent) :
    QDialog(parent),
    ui(new Ui::EquilateralTriangleParametersDialog)
{
    ui->setupUi(this);
}

EquilateralTriangleParametersDialog::~EquilateralTriangleParametersDialog()
{
    delete ui;
}

double EquilateralTriangleParametersDialog::GetLength() const
{
    return ui->lengthDoubleSpinBox->value();
}

QString EquilateralTriangleParametersDialog::GetPlane() const
{
    return ui->planeComboBox->currentText();
}
