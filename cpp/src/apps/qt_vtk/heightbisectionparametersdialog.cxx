#include "heightbisectionparametersdialog.h"
#include "ui_heightbisectionparametersdialog.h"

HeightBisectionParametersDialog::HeightBisectionParametersDialog(QWidget *parent) :
    QDialog(parent),
    ui(new Ui::HeightBisectionParametersDialog)
{
    ui->setupUi(this);
}

HeightBisectionParametersDialog::~HeightBisectionParametersDialog()
{
    delete ui;
}

double HeightBisectionParametersDialog::GetBottomHeight() const
{
    return ui->bottomHeightDoubleSpinBox->value();
}

bool HeightBisectionParametersDialog::GetBoundary() const
{
    return ui->boundaryCheckBox->isChecked();
}

double HeightBisectionParametersDialog::GetTopHeight() const
{
    return ui->topHeightDoubleSpinBox->value();
}
