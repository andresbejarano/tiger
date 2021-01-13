#ifndef EDGEDIRECTIONPARAMETERSDIALOG_H
#define EDGEDIRECTIONPARAMETERSDIALOG_H

#include <QDialog>

namespace Ui {
class EdgeDirectionParametersDialog;
}

class EdgeDirectionParametersDialog : public QDialog
{
    Q_OBJECT

public:

    /*
    Constructor of the class.
    */
    explicit EdgeDirectionParametersDialog(QWidget *parent = nullptr);

    /*
    Destructor of the class.
    */
    ~EdgeDirectionParametersDialog();

    /*
    Returns the selected initial direction value.
    @return QString The selected initial direction value.
    */
    QString GetInitialDirection() const;

private:

    Ui::EdgeDirectionParametersDialog *ui;
};

#endif // EDGEDIRECTIONPARAMETERSDIALOG_H
