#ifndef BENDEDSQUAREPARAMETERSDIALOG_H
#define BENDEDSQUAREPARAMETERSDIALOG_H

#include <QDialog>

namespace Ui {
class BendedSquareParametersDialog;
}

class BendedSquareParametersDialog : public QDialog
{
    Q_OBJECT

public:

    /*
    Constructor of the class.
    */
    explicit BendedSquareParametersDialog(QWidget *parent = nullptr);

    /*
    Destructor of the class.
    */
    ~BendedSquareParametersDialog();

    /*
    Returns the bend direction value from the dialog.
    @return QString The bend direction.
    */
    QString GetBendDirection() const;

    /*
    Returns the length value from the dialog.
    @return double The length value.
    */
    double GetLength() const;

    /*
    Returns the length segments value from the dialog.
    @return size_t The length segments value.
    */
    size_t GetLengthSegments() const;

    /*
    Returns the radial segments value from the dialog.
    @returns size_t The radial segments value.
    */
    size_t GetRadialSegments() const;

    /*
    Returns the radius value from the dialog.
    @return double The radius value.
    */
    double GetRadius() const;

private:

    Ui::BendedSquareParametersDialog *ui;
};

#endif // BENDEDSQUAREPARAMETERSDIALOG_H
