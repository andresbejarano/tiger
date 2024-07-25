#include "truncatepiecesparametersdialog.h"
#include "ui_truncatepiecesparametersdialog.h"

TruncatePiecesParametersDialog::TruncatePiecesParametersDialog(QWidget *parent) :
    QDialog(parent),
    ui(new Ui::TruncatePiecesParametersDialog)
{
    ui->setupUi(this);
}

TruncatePiecesParametersDialog::~TruncatePiecesParametersDialog()
{
    delete ui;
}

double TruncatePiecesParametersDialog::GetExtrados() const
{
    return ui->extradosDoubleSpinBox->value();
}

double TruncatePiecesParametersDialog::GetIntrados() const
{
    return ui->intradosDoubleSpinBox->value();
}
