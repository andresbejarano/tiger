#ifndef TILTINGANGLEPARAMETERSDIALOG_H
#define TILTINGANGLEPARAMETERSDIALOG_H

#include <QDialog>

namespace Ui {
class TiltingAngleParametersDialog;
}

class TiltingAngleParametersDialog : public QDialog
{
    Q_OBJECT

public:

    /*
    Constructor of the class.
    */
    explicit TiltingAngleParametersDialog(QWidget *parent = nullptr);

    /*
    Destructor of the class.
    */
    ~TiltingAngleParametersDialog();

    /*
    Returns the angle value for the Tilting Angle method.
    @return double The angle value for the Tilting Angle method.
    */
    double GetAngle() const;

    /*
    Returns the selected angle unit value.
    @return QString The selected angle unit value.
    */
    QString GetUnit() const;

private:

    Ui::TiltingAngleParametersDialog *ui;
};

#endif // TILTINGANGLEPARAMETERSDIALOG_H
