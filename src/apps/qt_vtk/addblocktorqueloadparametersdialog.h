#ifndef ADDBLOCKTORQUELOADPARAMETERSDIALOG_H
#define ADDBLOCKTORQUELOADPARAMETERSDIALOG_H

#include <QDialog>

namespace Ui 
{
    class AddBlockTorqueLoadParametersDialog;
}

class AddBlockTorqueLoadParametersDialog : public QDialog
{
    Q_OBJECT

public:

    // 
    // 
    // 
    explicit AddBlockTorqueLoadParametersDialog(QWidget *parent = nullptr);

    // 
    // 
    // 
    ~AddBlockTorqueLoadParametersDialog();

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
    Ui::AddBlockTorqueLoadParametersDialog *ui;
};

#endif // ADDBLOCKTORQUELOADPARAMETERSDIALOG_H
