#include "platonicsolidparametersdialog.h"
#include "ui_platonicsolidparametersdialog.h"

PlatonicSolidParametersDialog::PlatonicSolidParametersDialog(QWidget *parent) :
    QDialog(parent),
    ui(new Ui::PlatonicSolidParametersDialog)
{
    ui->setupUi(this);
}

PlatonicSolidParametersDialog::~PlatonicSolidParametersDialog()
{
    delete ui;
}

double PlatonicSolidParametersDialog::GetRadius() const
{
    return ui->radiusDoubleSpinBox->value();
}

QString PlatonicSolidParametersDialog::GetSolid() const
{
    return ui->solidComboBox->currentText();
}
