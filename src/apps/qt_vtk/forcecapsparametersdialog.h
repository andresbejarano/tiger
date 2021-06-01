#ifndef FORCECAPSPARAMETERSDIALOG_H
#define FORCECAPSPARAMETERSDIALOG_H

#include <QDialog>

namespace Ui 
{
    class ForceCapsParametersDialog;
}

class ForceCapsParametersDialog : public QDialog
{
    Q_OBJECT

public:

    //
    //
    //
    explicit ForceCapsParametersDialog(QWidget *parent = nullptr);

    //
    //
    //
    ~ForceCapsParametersDialog();

    //
    //
    //
    bool GetNoCaps() const;

    //
    //
    //
    double GetMinCap() const;

    //
    //
    //
    double GetMaxCap() const;

private:

    Ui::ForceCapsParametersDialog *ui;

private slots:

    void toggleCapsFields(int checkState);
};

#endif // FORCECAPSPARAMETERSDIALOG_H
