#ifndef TOGGLEBLOCKSPARAMETERSDIALOG_H
#define TOGGLEBLOCKSPARAMETERSDIALOG_H

#include <QDialog>

namespace Ui 
{
    class ToggleBlocksParametersDialog;
}

class ToggleBlocksParametersDialog : public QDialog
{
    Q_OBJECT

public:

    // 
    // Constructor of the class.
    // 
    explicit ToggleBlocksParametersDialog(QWidget *parent = nullptr);

    // 
    // Destructor of the class.
    // 
    ~ToggleBlocksParametersDialog();

    // 
    // 
    // 
    QString GetAction() const;

    // 
    // 
    // 
    QString GetBlockIndices() const;

    // 
    // 
    // 
    void SetBlockIndicesText(QString text);

private:

    Ui::ToggleBlocksParametersDialog *ui;
};

#endif // TOGGLEBLOCKSPARAMETERSDIALOG_H
