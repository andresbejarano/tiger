#ifndef REGULARPOLYGONPARAMETERSDIALOG_H
#define REGULARPOLYGONPARAMETERSDIALOG_H

#include <QDialog>

namespace Ui {
class RegularPolygonParametersDialog;
}

class RegularPolygonParametersDialog : public QDialog
{
    Q_OBJECT

public:

    /*
    Constructor of the class.
    */
    explicit RegularPolygonParametersDialog(QWidget *parent = nullptr);

    /*
    Destructor of the class.
    */
    ~RegularPolygonParametersDialog();

    /*
    Returns the length of the sides of the regular polygon.
    @return double The length of the sides of the regular polygon.
    */
    double GetLength() const;

    /*
    Returns the selected plane value.
    @return QString The selected plane value.
    */
    QString GetPlane() const;

    /*
    Returns the number of sides of the regular polygon.
    @return size_t The number of sides of the regular polygon.
    */
    size_t GetSides() const;

private:

    Ui::RegularPolygonParametersDialog *ui;
};

#endif // REGULARPOLYGONPARAMETERSDIALOG_H
