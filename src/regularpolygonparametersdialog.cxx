#include "regularpolygonparametersdialog.h"
#include "ui_regularpolygonparametersdialog.h"

RegularPolygonParametersDialog::RegularPolygonParametersDialog(QWidget *parent) :
    QDialog(parent),
    ui(new Ui::RegularPolygonParametersDialog)
{
    ui->setupUi(this);
}

RegularPolygonParametersDialog::~RegularPolygonParametersDialog()
{
    delete ui;
}

double RegularPolygonParametersDialog::GetLength() const
{
    return ui->lengthDoubleSpinBox->value();
}

QString RegularPolygonParametersDialog::GetPlane() const
{
    return ui->planeComboBox->currentText();
}

size_t RegularPolygonParametersDialog::GetSides() const
{
    return ui->sidesSpinBox->value();
}
