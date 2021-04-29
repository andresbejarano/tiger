#include "normalizetessellationparametersdialog.h"
#include "ui_normalizetessellationparametersdialog.h"

NormalizeTessellationParametersDialog::NormalizeTessellationParametersDialog(QWidget *parent) :
    QDialog(parent),
    ui(new Ui::NormalizeTessellationParametersDialog)
{
    ui->setupUi(this);
}

NormalizeTessellationParametersDialog::~NormalizeTessellationParametersDialog()
{
    delete ui;
}

double NormalizeTessellationParametersDialog::GetRadius() const
{
    return ui->radiusDoubleSpinBox->value();
}
