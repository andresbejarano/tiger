#include "archimedeansolidparametersdialog.h"
#include "ui_archimedeansolidparametersdialog.h"

ArchimedeanSolidParametersDialog::ArchimedeanSolidParametersDialog(QWidget *parent) :
    QDialog(parent),
    ui(new Ui::ArchimedeanSolidParametersDialog)
{
    ui->setupUi(this);
}

ArchimedeanSolidParametersDialog::~ArchimedeanSolidParametersDialog()
{
    delete ui;
}

QString ArchimedeanSolidParametersDialog::GetSolid() const
{
    return ui->solidComboBox->currentText();
}

double ArchimedeanSolidParametersDialog::GetRadius() const
{
    return ui->radiusDoubleSpinBox->value();
}
