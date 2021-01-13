#include <qmessagebox.h>
#include "squareparametersdialog.h"
#include "ui_squareparametersdialog.h"

SquareParametersDialog::SquareParametersDialog(QWidget *parent) :
    QDialog(parent),
    ui(new Ui::SquareParametersDialog)
{
    ui->setupUi(this);

	// Verify the correctness of the dialog values when user clicks OK
	connect(ui->buttonBox, &QDialogButtonBox::accepted, this, &SquareParametersDialog::verify);
}

SquareParametersDialog::~SquareParametersDialog()
{
    delete ui;
}

double SquareParametersDialog::GetLength() const
{
	return ui->lengthDoubleSpinBox->value();
}

QString SquareParametersDialog::GetPlane() const
{
	return ui->planeComboBox->currentText();
}

void SquareParametersDialog::verify() 
{
	// 
	if (ui->lengthDoubleSpinBox->value() > 0.0)
	{
		accept();
		return;
	}

	// 
	QMessageBox::warning(this, "Invalid Value", "Length value must be nonnegative.");
}