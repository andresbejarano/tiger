#ifndef SQUAREPARAMETERSDIALOG_H
#define SQUAREPARAMETERSDIALOG_H

#include <QDialog>

namespace Ui {
class SquareParametersDialog;
}

class SquareParametersDialog : public QDialog
{
    Q_OBJECT

public:

	/*
	Constructor of the class.
	*/
    explicit SquareParametersDialog(QWidget *parent = nullptr);

	/*
	Destructor of the class.
	*/
    ~SquareParametersDialog();

	/*
	Returns the length value of the dialog.
	@return double The length value of the dialog.
	*/
	double GetLength() const;

	/*
	Returns the plane value of the dialog.
	@return QString The plane value of the dialog.
	*/
	QString GetPlane() const;

public slots:

	/*
	Verifies the correctness of the input values of the dialog.
	*/
	void verify();

private:
    Ui::SquareParametersDialog *ui;
};

#endif // SQUAREPARAMETERSDIALOG_H
