#include "forcecapsparametersdialog.h"
#include "ui_forcecapsparametersdialog.h"

ForceCapsParametersDialog::ForceCapsParametersDialog(QWidget *parent) :
    QDialog(parent),
    ui(new Ui::ForceCapsParametersDialog)
{
    ui->setupUi(this);

    // Connect the no caps checkbox state signal with the function that toggles
    // the min and max force caps fields
    QObject::connect(
        ui->noCapsCheckBox, 
        SIGNAL(stateChanged(int)),
        this, 
        SLOT(toggleCapsFields(int)));
}

ForceCapsParametersDialog::~ForceCapsParametersDialog()
{
    delete ui;
}

bool ForceCapsParametersDialog::GetNoCaps() const
{
    return ui->noCapsCheckBox->isChecked();
}

double ForceCapsParametersDialog::GetMinCap() const
{
    return ui->minCapDoubleSpinBox->value();
}

double ForceCapsParametersDialog::GetMaxCap() const
{
    return ui->maxCapDoubleSpinBox->value();
}

void ForceCapsParametersDialog::toggleCapsFields(int checkState)
{
    bool useCapsChecked = (checkState == Qt::Unchecked);

    ui->minCapDoubleSpinBox->setEnabled(useCapsChecked);
    ui->maxCapDoubleSpinBox->setEnabled(useCapsChecked);
}
