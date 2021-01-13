#ifndef PARABOLOIDPARAMETERSDIALOG_H
#define PARABOLOIDPARAMETERSDIALOG_H

#include <QDialog>

namespace Ui {
class ParaboloidParametersDialog;
}

class ParaboloidParametersDialog : public QDialog
{
    Q_OBJECT

public:

    explicit ParaboloidParametersDialog(QWidget *parent = nullptr);

    ~ParaboloidParametersDialog();

    double GetWidth() const;

    double GetHeight() const;

    size_t GetWidthSegments() const;

    size_t GetHeightSegments() const;

    QString GetPlane() const;

    double GetA() const;

    double GetB() const;

    QString GetType() const;

private:
    Ui::ParaboloidParametersDialog *ui;
};

#endif // PARABOLOIDPARAMETERSDIALOG_H
