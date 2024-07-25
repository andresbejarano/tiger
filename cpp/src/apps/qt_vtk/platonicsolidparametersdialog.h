#ifndef PLATONICSOLIDPARAMETERSDIALOG_H
#define PLATONICSOLIDPARAMETERSDIALOG_H

#include <QDialog>

namespace Ui {
class PlatonicSolidParametersDialog;
}

class PlatonicSolidParametersDialog : public QDialog
{
    Q_OBJECT

public:

    /*
    Constructor of the class.
    */
    explicit PlatonicSolidParametersDialog(QWidget *parent = nullptr);

    /*
    Destructor of the class.
    */
    ~PlatonicSolidParametersDialog();

    /*
    Returns the radius of the solid.
    @return double The radius of the solid.
    */
    double GetRadius() const;

    /*
    Returns the selected Platonic solid.
    @return QString The selected Platonic solid.
    */
    QString GetSolid() const;

private:

    Ui::PlatonicSolidParametersDialog *ui;
};

#endif // PLATONICSOLIDPARAMETERSDIALOG_H
