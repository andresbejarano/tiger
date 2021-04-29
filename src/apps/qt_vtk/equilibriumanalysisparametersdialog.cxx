#include "equilibriumanalysisparametersdialog.h"
#include "ui_equilibriumanalysisparametersdialog.h"

EquilibriumAnalysisParametersDialog::EquilibriumAnalysisParametersDialog(QWidget *parent) :
    QDialog(parent),
    ui(new Ui::EquilibriumAnalysisParametersDialog)
{
    ui->setupUi(this);
}

EquilibriumAnalysisParametersDialog::~EquilibriumAnalysisParametersDialog()
{
    delete ui;
}

double EquilibriumAnalysisParametersDialog::GetDensity() const
{
    return ui->densityDoubleSpinBox->value();
}

double EquilibriumAnalysisParametersDialog::GetFrictionCoefficient() const
{
    return ui->frictionCoefficientDoubleSpinBox->value();
}

QString EquilibriumAnalysisParametersDialog::GetLengthUnit() const
{
    return ui->lengthUnitComboBox->currentText();
}

double EquilibriumAnalysisParametersDialog::GetCompressionWeight() const
{
    return ui->compressionWeightDoubleSpinBox->value();
}

double EquilibriumAnalysisParametersDialog::GetTensionWeight() const
{
    return ui->tensionWeightDoubleSpinBox->value();
}

double EquilibriumAnalysisParametersDialog::GetUTangentialWeight() const
{
    return ui->uTangentialWeightDoubleSpinBox->value();
}

double EquilibriumAnalysisParametersDialog::GetVTangentialWeight() const
{
    return ui->vTangentialWeightDoubleSpinBox->value();
}

bool EquilibriumAnalysisParametersDialog::GetNormalize() const
{
    return ui->normalizeForcesCheckBox->isChecked();
}

bool EquilibriumAnalysisParametersDialog::GetVerbose() const
{
    return ui->verboseCheckBox->isChecked();
}

bool EquilibriumAnalysisParametersDialog::GetFiles() const
{
    return ui->filesCheckBox->isChecked();
}
