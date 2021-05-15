#ifndef RESETBLOCKLOADSPARAMETERSDIALOG_H
#define RESETBLOCKLOADSPARAMETERSDIALOG_H

#include <QDialog>

namespace Ui 
{
    class ResetBlockLoadsParametersDialog;
}

class ResetBlockLoadsParametersDialog : public QDialog
{
    Q_OBJECT

public:

    // 
    // 
    // 
    explicit ResetBlockLoadsParametersDialog(QWidget *parent = nullptr);

    // 
    // 
    // 
    ~ResetBlockLoadsParametersDialog();

    // 
    // 
    // 
    QString GetBlockIndices() const;

    // 
    // 
    // 
    QString GetLoadType() const;

private:

    Ui::ResetBlockLoadsParametersDialog *ui;
};

#endif // RESETBLOCKLOADSPARAMETERSDIALOG_H
