#include "rectangleparametersdialog.h"
#include "ui_rectangleparametersdialog.h"

RectangleParametersDialog::RectangleParametersDialog(QWidget *parent) :
    QDialog(parent),
    ui(new Ui::RectangleParametersDialog)
{
    ui->setupUi(this);
}

RectangleParametersDialog::~RectangleParametersDialog()
{
    delete ui;
}

double RectangleParametersDialog::GetHeight() const
{
    return ui->heightDoubleSpinBox->value();
}

QString RectangleParametersDialog::GetPlane() const
{
    return ui->planeComboBox->currentText();
}

double RectangleParametersDialog::GetWidth() const
{
    return ui->widthDoubleSpinBox->value();
}
