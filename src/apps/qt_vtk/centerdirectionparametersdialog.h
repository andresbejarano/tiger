#ifndef CENTERDIRECTIONPARAMETERSDIALOG_H
#define CENTERDIRECTIONPARAMETERSDIALOG_H

#include <QDialog>

namespace Ui {
class CenterDirectionParametersDialog;
}

class CenterDirectionParametersDialog : public QDialog
{
    Q_OBJECT

public:

    /*
    Constructor of the class.
    */
    explicit CenterDirectionParametersDialog(QWidget *parent = nullptr);

    /*
    Destructor of the class.
    */
    ~CenterDirectionParametersDialog();

    /*
    */
    QString GetCenter() const;

    /*
    */
    QString GetDirection() const;

private:

    Ui::CenterDirectionParametersDialog *ui;
};

#endif // CENTERDIRECTIONPARAMETERSDIALOG_H
