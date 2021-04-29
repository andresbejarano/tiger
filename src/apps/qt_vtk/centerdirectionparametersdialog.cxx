#include "centerdirectionparametersdialog.h"
#include "ui_centerdirectionparametersdialog.h"

CenterDirectionParametersDialog::CenterDirectionParametersDialog(QWidget *parent) :
    QDialog(parent),
    ui(new Ui::CenterDirectionParametersDialog)
{
    ui->setupUi(this);
}

CenterDirectionParametersDialog::~CenterDirectionParametersDialog()
{
    delete ui;
}

QString CenterDirectionParametersDialog::GetCenter() const
{
    return ui->centerComboBox->currentText();
}

QString CenterDirectionParametersDialog::GetDirection() const
{
    return ui->directionComboBox->currentText();
}
