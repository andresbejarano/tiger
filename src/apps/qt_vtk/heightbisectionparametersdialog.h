#ifndef HEIGHTBISECTIONPARAMETERSDIALOG_H
#define HEIGHTBISECTIONPARAMETERSDIALOG_H

#include <QDialog>

namespace Ui {
class HeightBisectionParametersDialog;
}

class HeightBisectionParametersDialog : public QDialog
{
    Q_OBJECT

public:

    /*
    Constructor of the class.
    */
    explicit HeightBisectionParametersDialog(QWidget *parent = nullptr);

    /*
    Destructor of the class.
    */
    ~HeightBisectionParametersDialog();

    /*
    Returns the bottom height value for the Height-Bisection method.
    @return double The bottom height value value for the Height-Bisection method.
    */
    double GetBottomHeight() const;

    /*
    Returns the indicator for generating the pieces at the boundary.
    @return bool The indicator for generating the pieces at the boundary.
    */
    bool GetBoundary() const;

    /*
    Returns the top height value for the Height-Bisection method.
    @return double The top height value for the Height-Bisection method.
    */
    double GetTopHeight() const;

private:

    Ui::HeightBisectionParametersDialog *ui;
};

#endif // HEIGHTBISECTIONPARAMETERSDIALOG_H
