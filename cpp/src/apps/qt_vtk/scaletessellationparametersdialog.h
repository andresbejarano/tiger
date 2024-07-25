#ifndef SCALETESSELLATIONPARAMETERSDIALOG_H
#define SCALETESSELLATIONPARAMETERSDIALOG_H

#include <QDialog>

namespace Ui {
class ScaleTessellationParametersDialog;
}

class ScaleTessellationParametersDialog : public QDialog
{
    Q_OBJECT

public:
    explicit ScaleTessellationParametersDialog(QWidget *parent = nullptr);
    ~ScaleTessellationParametersDialog();

    double GetFactor() const;

private:
    Ui::ScaleTessellationParametersDialog *ui;
};

#endif // SCALETESSELLATIONPARAMETERSDIALOG_H
