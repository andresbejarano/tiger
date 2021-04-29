#include "rotatetessellationparametersdialog.h"
#include "ui_rotatetessellationparametersdialog.h"

RotateTessellationParametersDialog::RotateTessellationParametersDialog(QWidget *parent) :
    QDialog(parent),
    ui(new Ui::RotateTessellationParametersDialog)
{
    ui->setupUi(this);
}

RotateTessellationParametersDialog::~RotateTessellationParametersDialog()
{
    delete ui;
}

double RotateTessellationParametersDialog::GetX() const
{
    return ui->xDoubleSpinBox->value();
}

double RotateTessellationParametersDialog::GetY() const
{
    return ui->yDoubleSpinBox->value();
}

double RotateTessellationParametersDialog::GetZ() const
{
    return ui->zDoubleSpinBox->value();
}

double RotateTessellationParametersDialog::GetAngle() const
{
    return ui->angleDoubleSpinBox->value();
}
