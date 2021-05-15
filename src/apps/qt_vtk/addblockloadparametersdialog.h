#ifndef ADDBLOCKLOADPARAMETERSDIALOG_H
#define ADDBLOCKLOADPARAMETERSDIALOG_H

#include <QDialog>

namespace Ui 
{
    class AddBlockLoadParametersDialog;
}

class AddBlockLoadParametersDialog : public QDialog
{
    Q_OBJECT

public:

    // 
    // 
    // 
    explicit AddBlockLoadParametersDialog(QWidget *parent = nullptr);

    // 
    // 
    // 
    ~AddBlockLoadParametersDialog();

    // 
    // 
    // 
    QString GetBlockIndices() const;

    // 
    // 
    // 
    QString GetLoadType() const;

    // 
    // 
    // 
    double GetX() const;

    // 
    // 
    // 
    double GetY() const;

    // 
    // 
    // 
    double GetZ() const;

private:

    Ui::AddBlockLoadParametersDialog *ui;
};

#endif // ADDBLOCKLOADPARAMETERSDIALOG_H
