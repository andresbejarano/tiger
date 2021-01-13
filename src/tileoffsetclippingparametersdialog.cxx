#include "tileoffsetclippingparametersdialog.h"
#include "ui_tileoffsetclippingparametersdialog.h"

TileOffsetClippingParametersDialog::TileOffsetClippingParametersDialog(QWidget *parent) :
    QDialog(parent),
    ui(new Ui::TileOffsetClippingParametersDialog)
{
    ui->setupUi(this);
}

TileOffsetClippingParametersDialog::~TileOffsetClippingParametersDialog()
{
    delete ui;
}

double TileOffsetClippingParametersDialog::GetExtrados() const
{
    return ui->extradosDoubleSpinBox->value();
}

double TileOffsetClippingParametersDialog::GetIntrados() const
{
    return ui->intradosDoubleSpinBox->value();
}
