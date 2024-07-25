#include "toggleblocksparametersdialog.h"
#include "ui_toggleblocksparametersdialog.h"

ToggleBlocksParametersDialog::ToggleBlocksParametersDialog(QWidget *parent) :
    QDialog(parent),
    ui(new Ui::ToggleBlocksParametersDialog)
{
    ui->setupUi(this);
}

ToggleBlocksParametersDialog::~ToggleBlocksParametersDialog()
{
    delete ui;
}

QString ToggleBlocksParametersDialog::GetAction() const
{
    return ui->actionComboBox->currentText();
}

QString ToggleBlocksParametersDialog::GetBlockIndices() const
{
    return ui->blockIndicesLineEdit->text();
}

void ToggleBlocksParametersDialog::SetBlockIndicesText(QString text)
{
    ui->blockIndicesLineEdit->setText(text);
}
