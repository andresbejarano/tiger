#ifndef CYLINDERPARAMETERSDIALOG_H
#define CYLINDERPARAMETERSDIALOG_H

#include <QDialog>

namespace Ui {
class CylinderParametersDialog;
}

class CylinderParametersDialog : public QDialog
{
    Q_OBJECT

public:

    /*
    Constructor of the class.
    */
    explicit CylinderParametersDialog(QWidget *parent = nullptr);

    /*
    Destructor of the class.
    */
    ~CylinderParametersDialog();

    /*
    Returns the selected axis value.
    @return QString The selected axis value.
    */
    QString GetAxis() const;

    /*
    Returns the length of the cylinder.
    @return double The length of the cylinder.
    */
    double GetLength() const;

    /*
    Returns the number of length segments of the cylinder.
    @return size_t The number of length segments of the cylinder.
    */
    size_t GetLengthSegments() const;

    /*
    Returns the number of radial segments of the cylinder.
    @return size_t The number of radial segments of the cylinder.
    */
    size_t GetRadialSegments() const;

    /*
    Returns the radius of the cylinder.
    @return double The radius of the cylinder.
    */
    double GetRadius() const;

private:

    Ui::CylinderParametersDialog *ui;
};

#endif // CYLINDERPARAMETERSDIALOG_H
