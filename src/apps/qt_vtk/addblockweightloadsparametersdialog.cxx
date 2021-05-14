#include "addblockweightloadsparametersdialog.h"
#include "ui_addblockweightloadsparametersdialog.h"

AddBlockWeightLoadsParametersDialog::AddBlockWeightLoadsParametersDialog(QWidget *parent) :
    QDialog(parent),
    ui(new Ui::AddBlockWeightLoadsParametersDialog)
{
    ui->setupUi(this);
}

AddBlockWeightLoadsParametersDialog::~AddBlockWeightLoadsParametersDialog()
{
    delete ui;
}

QString AddBlockWeightLoadsParametersDialog::GetBlockIndices() const
{
    return ui->blockIndicesLineEdit->text();
}

double AddBlockWeightLoadsParametersDialog::GetDensity() const
{
    return ui->blockDensityDoubleSpinBox->value();
}

QString AddBlockWeightLoadsParametersDialog::GetUnits() const
{
    return ui->blockUnitComboBox->currentText();
}
