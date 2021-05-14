#include "addblocktorqueloadparametersdialog.h"
#include "ui_addblocktorqueloadparametersdialog.h"

AddBlockTorqueLoadParametersDialog::AddBlockTorqueLoadParametersDialog(QWidget *parent) :
    QDialog(parent),
    ui(new Ui::AddBlockTorqueLoadParametersDialog)
{
    ui->setupUi(this);
}

AddBlockTorqueLoadParametersDialog::~AddBlockTorqueLoadParametersDialog()
{
    delete ui;
}

QString AddBlockTorqueLoadParametersDialog::GetBlockIndices() const
{
    return ui->blockIndicesLineEdit->text();
}

QString AddBlockTorqueLoadParametersDialog::GetLoadType() const 
{
    return ui->loadTypeComboBox->currentText();
}

double AddBlockTorqueLoadParametersDialog::GetX() const
{
    return ui->xDoubleSpinBox->value();
}

double AddBlockTorqueLoadParametersDialog::GetY() const
{
    return ui->yDoubleSpinBox->value();
}

double AddBlockTorqueLoadParametersDialog::GetZ() const
{
    return ui->zDoubleSpinBox->value();
}
