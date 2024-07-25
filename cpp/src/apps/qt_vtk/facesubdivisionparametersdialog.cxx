#include "facesubdivisionparametersdialog.h"
#include "ui_facesubdivisionparametersdialog.h"

FaceSubdivisionParametersDialog::FaceSubdivisionParametersDialog(QWidget *parent) :
    QDialog(parent),
    ui(new Ui::FaceSubdivisionParametersDialog)
{
    ui->setupUi(this);
}

FaceSubdivisionParametersDialog::~FaceSubdivisionParametersDialog()
{
    delete ui;
}

QString FaceSubdivisionParametersDialog::GetCenter() const
{
    return ui->centerComboBox->currentText();
}

QString FaceSubdivisionParametersDialog::GetType() const
{
    return ui->typeComboBox->currentText();
}
