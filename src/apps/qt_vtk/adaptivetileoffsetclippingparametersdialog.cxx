#include "adaptivetileoffsetclippingparametersdialog.h"
#include "ui_adaptivetileoffsetclippingparametersdialog.h"

AdaptiveTileOffsetClippingParametersDialog::AdaptiveTileOffsetClippingParametersDialog(QWidget *parent) :
    QDialog(parent),
    ui(new Ui::AdaptiveTileOffsetClippingParametersDialog)
{
    ui->setupUi(this);
}

AdaptiveTileOffsetClippingParametersDialog::~AdaptiveTileOffsetClippingParametersDialog()
{
    delete ui;
}

QString AdaptiveTileOffsetClippingParametersDialog::GetBottomFunction() const
{
    return ui->bottomFunctionComboBox->currentText();
}

QString AdaptiveTileOffsetClippingParametersDialog::GetTopFunction() const
{
    return ui->topFunctionComboBox->currentText();
}
