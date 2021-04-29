#ifndef NORMALIZETESSELLATIONPARAMETERSDIALOG_H
#define NORMALIZETESSELLATIONPARAMETERSDIALOG_H

#include <QDialog>

namespace Ui {
class NormalizeTessellationParametersDialog;
}

class NormalizeTessellationParametersDialog : public QDialog
{
    Q_OBJECT

public:
    explicit NormalizeTessellationParametersDialog(QWidget *parent = nullptr);
    ~NormalizeTessellationParametersDialog();

    double GetRadius() const;

private:
    Ui::NormalizeTessellationParametersDialog *ui;
};

#endif // NORMALIZETESSELLATIONPARAMETERSDIALOG_H
