#include "scaletessellationparametersdialog.h"
#include "ui_scaletessellationparametersdialog.h"

ScaleTessellationParametersDialog::ScaleTessellationParametersDialog(QWidget *parent) :
    QDialog(parent),
    ui(new Ui::ScaleTessellationParametersDialog)
{
    ui->setupUi(this);
}

ScaleTessellationParametersDialog::~ScaleTessellationParametersDialog()
{
    delete ui;
}

double ScaleTessellationParametersDialog::GetFactor() const
{
    return ui->factorDoubleSpinBox->value();
}
