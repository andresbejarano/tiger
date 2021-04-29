#ifndef TRUNCATEPIECESPARAMETERSDIALOG_H
#define TRUNCATEPIECESPARAMETERSDIALOG_H

#include <QDialog>

namespace Ui {
class TruncatePiecesParametersDialog;
}

class TruncatePiecesParametersDialog : public QDialog
{
    Q_OBJECT

public:

    /*
    Constructor of the class.
    */
    explicit TruncatePiecesParametersDialog(QWidget *parent = nullptr);

    /*
    Destructor of the class.
    */
    ~TruncatePiecesParametersDialog();

    /*
    Returns the extrados value of the dialog.
    @return double The extrados value of the dialog.
    */
    double GetExtrados() const;

    /*
    Returns the intrados value of the dialog.
    @return double The intrados value of the dialog.
    */
    double GetIntrados() const;

private:
    Ui::TruncatePiecesParametersDialog *ui;
};

#endif // TRUNCATEPIECESPARAMETERSDIALOG_H
