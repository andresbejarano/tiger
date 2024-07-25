#include "addblockloadparametersdialog.h"
#include "ui_addblockloadparametersdialog.h"

AddBlockLoadParametersDialog::AddBlockLoadParametersDialog(QWidget *parent) :
    QDialog(parent),
    ui(new Ui::AddBlockLoadParametersDialog)
{
    ui->setupUi(this);
}

AddBlockLoadParametersDialog::~AddBlockLoadParametersDialog()
{
    delete ui;
}

QString AddBlockLoadParametersDialog::GetBlockIndices() const
{
    return ui->blockIndicesLineEdit->text();
}

QString AddBlockLoadParametersDialog::GetLoadType() const
{
    return ui->loadTypeComboBox->currentText();
}

double AddBlockLoadParametersDialog::GetX() const
{
    return ui->xDoubleSpinBox->value();
}

double AddBlockLoadParametersDialog::GetY() const
{
    return ui->yDoubleSpinBox->value();
}

double AddBlockLoadParametersDialog::GetZ() const
{
    return ui->zDoubleSpinBox->value();
}
