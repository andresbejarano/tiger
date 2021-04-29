#ifndef CONEPARAMETERSDIALOG_H
#define CONEPARAMETERSDIALOG_H

#include <QDialog>

namespace Ui {
class ConeParametersDialog;
}

class ConeParametersDialog : public QDialog
{
    Q_OBJECT

public:

    /*
    Constructor of the class.
    */
    explicit ConeParametersDialog(QWidget *parent = nullptr);

    /*
    Destructor of the class.
    */
    ~ConeParametersDialog();

    /*
    Returns the name of the selected axis value.
    @return QString The name of the selected axis value.
    */
    QString GetAxis() const;

    /*
    Returns the length of the cone.
    @return double The length of the cone.
    */
    double GetLength() const;

    /*
    Returns the number of length segments of the cone.
    @return size_t The number of length segments of the cone.
    */
    size_t GetLengthSegments() const;

    /*
    Returns the number of radial segments of the cone.
    @return size_t The number of radial segments of the cone.
    */
    size_t GetRadialSegments() const;

    /*
    Returns the radius of the cone.
    @return double The radius of the cone.
    */
    double GetRadius() const;

private:

    Ui::ConeParametersDialog *ui;
};

#endif // CONEPARAMETERSDIALOG_H
