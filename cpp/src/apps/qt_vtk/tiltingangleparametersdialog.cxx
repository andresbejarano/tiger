#include "tiltingangleparametersdialog.h"
#include "ui_tiltingangleparametersdialog.h"

TiltingAngleParametersDialog::TiltingAngleParametersDialog(QWidget *parent) :
    QDialog(parent),
    ui(new Ui::TiltingAngleParametersDialog)
{
    ui->setupUi(this);
}

TiltingAngleParametersDialog::~TiltingAngleParametersDialog()
{
    delete ui;
}

double TiltingAngleParametersDialog::GetAngle() const
{
    return ui->angleDoubleSpinBox->value();
}

QString TiltingAngleParametersDialog::GetUnit() const
{
    return ui->unitComboBox->currentText();
}
