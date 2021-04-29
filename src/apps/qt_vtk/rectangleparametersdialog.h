#ifndef RECTANGLEPARAMETERSDIALOG_H
#define RECTANGLEPARAMETERSDIALOG_H

#include <QDialog>

namespace Ui {
class RectangleParametersDialog;
}

class RectangleParametersDialog : public QDialog
{
    Q_OBJECT

public:

    /*
    Constructor of the class.
    */
    explicit RectangleParametersDialog(QWidget *parent = nullptr);

    /*
    Destructor of the class.
    */
    ~RectangleParametersDialog();

    /*
    Returns the height of the rectangle.
    @return double The height of the rectangle.
    */
    double GetHeight() const;

    /*
    Returns the selected plane value.
    @returns QString The selected plane value.
    */
    QString GetPlane() const;

    /*
    Returns the width of the rectangle.
    @return double The width of the rectangle.
    */
    double GetWidth() const;

private:

    Ui::RectangleParametersDialog *ui;
};

#endif // RECTANGLEPARAMETERSDIALOG_H
