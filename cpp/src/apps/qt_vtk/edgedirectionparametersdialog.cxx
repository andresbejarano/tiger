#include "edgedirectionparametersdialog.h"
#include "ui_edgedirectionparametersdialog.h"

EdgeDirectionParametersDialog::EdgeDirectionParametersDialog(QWidget *parent) :
    QDialog(parent),
    ui(new Ui::EdgeDirectionParametersDialog)
{
    ui->setupUi(this);
}

EdgeDirectionParametersDialog::~EdgeDirectionParametersDialog()
{
    delete ui;
}

QString EdgeDirectionParametersDialog::GetInitialDirection() const
{
    return ui->initialDirectionComboBox->currentText();
}
