#ifndef EQUILATERALTRIANGLEPARAMETERSDIALOG_H
#define EQUILATERALTRIANGLEPARAMETERSDIALOG_H

#include <QDialog>

namespace Ui {
class EquilateralTriangleParametersDialog;
}

class EquilateralTriangleParametersDialog : public QDialog
{
    Q_OBJECT

public:

    /*
    Constructor of the class.
    */
    explicit EquilateralTriangleParametersDialog(QWidget *parent = nullptr);

    /*
    Destructor of the class.
    */
    ~EquilateralTriangleParametersDialog();

    /*
    Returns the length of the sides of the equilateral triangle.
    @return double The length of the sides of the equilateral triangle.
    */
    double GetLength() const;

    /*
    Returns the selected plane value.
    @return QString The selected plane value.
    */
    QString GetPlane() const;

private:

    Ui::EquilateralTriangleParametersDialog *ui;
};

#endif // EQUILATERALTRIANGLEPARAMETERSDIALOG_H
