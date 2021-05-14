#ifndef ADDBLOCKWEIGHTLOADSPARAMETERSDIALOG_H
#define ADDBLOCKWEIGHTLOADSPARAMETERSDIALOG_H

#include <QDialog>

namespace Ui 
{
    class AddBlockWeightLoadsParametersDialog;
}

class AddBlockWeightLoadsParametersDialog : public QDialog
{
    Q_OBJECT

public:

    // 
    //
    // 
    explicit AddBlockWeightLoadsParametersDialog(QWidget *parent = nullptr);

    // 
    // 
    // 
    ~AddBlockWeightLoadsParametersDialog();

    // 
    // 
    // 
    QString GetBlockIndices() const;

    // 
    // 
    // 
    double GetDensity() const;

    // 
    // 
    // 
    QString GetUnits() const;

private:
    Ui::AddBlockWeightLoadsParametersDialog *ui;
};

#endif // ADDBLOCKWEIGHTLOADSPARAMETERSDIALOG_H
