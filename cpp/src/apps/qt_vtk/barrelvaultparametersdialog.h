#ifndef BARRELVAULTPARAMETERSDIALOG_H
#define BARRELVAULTPARAMETERSDIALOG_H

#include <QDialog>

namespace Ui 
{
    class BarrelVaultParametersDialog;
}

class BarrelVaultParametersDialog : public QDialog
{
    Q_OBJECT

public:

    //
    // Constructor of the class.
    //
    explicit BarrelVaultParametersDialog(QWidget *parent = nullptr);

    //
    //
    //
    ~BarrelVaultParametersDialog();

    // 
    // 
    // 
    double GetLength() const;

    // 
    // 
    // 
    double GetRadius() const;

    // 
    // 
    // 
    size_t GetLengthSegments() const;

    // 
    // 
    // 
    size_t GetRadialSegments() const;

    // 
    // 
    // 
    QString GetAxis() const;

private:

    // 
    Ui::BarrelVaultParametersDialog *ui;
};

#endif // BARRELVAULTPARAMETERSDIALOG_H
