#ifndef ARCHIMEDEANSOLIDPARAMETERSDIALOG_H
#define ARCHIMEDEANSOLIDPARAMETERSDIALOG_H

#include <QDialog>

namespace Ui {
class ArchimedeanSolidParametersDialog;
}

class ArchimedeanSolidParametersDialog : public QDialog
{
    Q_OBJECT

public:

    /*
    Constructor of the class.
    */
    explicit ArchimedeanSolidParametersDialog(QWidget *parent = nullptr);

    /*
    Destructor of the class.
    */
    ~ArchimedeanSolidParametersDialog();

    /*
    Returns the name of the selected Archimedean solid.
    @return QString The name of the selected Archimedean solid.
    */
    QString GetSolid() const;

    /*
    Returns the radius of the solid.
    @return double The radius of the solid.
    */
    double GetRadius() const;

private:

    Ui::ArchimedeanSolidParametersDialog *ui;
};

#endif // ARCHIMEDEANSOLIDPARAMETERSDIALOG_H
