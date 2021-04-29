#ifndef ROTATETESSELLATIONPARAMETERSDIALOG_H
#define ROTATETESSELLATIONPARAMETERSDIALOG_H

#include <QDialog>

namespace Ui {
class RotateTessellationParametersDialog;
}

class RotateTessellationParametersDialog : public QDialog
{
    Q_OBJECT

public:
    explicit RotateTessellationParametersDialog(QWidget *parent = nullptr);
    ~RotateTessellationParametersDialog();

    double GetX() const;

    double GetY() const;

    double GetZ() const;

    double GetAngle() const;

private:
    Ui::RotateTessellationParametersDialog *ui;
};

#endif // ROTATETESSELLATIONPARAMETERSDIALOG_H
