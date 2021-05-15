#include "resetblockloadsparametersdialog.h"
#include "ui_resetblockloadsparametersdialog.h"

ResetBlockLoadsParametersDialog::ResetBlockLoadsParametersDialog(QWidget *parent) :
    QDialog(parent),
    ui(new Ui::ResetBlockLoadsParametersDialog)
{
    ui->setupUi(this);
}

ResetBlockLoadsParametersDialog::~ResetBlockLoadsParametersDialog()
{
    delete ui;
}

QString ResetBlockLoadsParametersDialog::GetBlockIndices() const
{
    return ui->blockIndicesLineEdit->text();
}

QString ResetBlockLoadsParametersDialog::GetLoadType() const
{
    return ui->resetLoadComboBox->currentText();
}

