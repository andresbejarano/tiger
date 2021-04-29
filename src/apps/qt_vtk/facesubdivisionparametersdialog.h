#ifndef FACESUBDIVISIONPARAMETERSDIALOG_H
#define FACESUBDIVISIONPARAMETERSDIALOG_H

#include <QDialog>

namespace Ui {
class FaceSubdivisionParametersDialog;
}

class FaceSubdivisionParametersDialog : public QDialog
{
    Q_OBJECT

public:

    /*
    Constructor of the class.
    */
    explicit FaceSubdivisionParametersDialog(QWidget *parent = nullptr);

    /*
    Destructor of the class.
    */
    ~FaceSubdivisionParametersDialog();

    /*
    Returns the selected center value.
    @return QString The selected center value.
    */
    QString GetCenter() const;

    /*
    Returns the selected type value.
    @return QString The selected type value.
    */
    QString GetType() const;

private:

    Ui::FaceSubdivisionParametersDialog *ui;
};

#endif // FACESUBDIVISIONPARAMETERSDIALOG_H
