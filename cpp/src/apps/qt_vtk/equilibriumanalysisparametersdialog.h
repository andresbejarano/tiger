#ifndef EQUILIBRIUMANALYSISPARAMETERSDIALOG_H
#define EQUILIBRIUMANALYSISPARAMETERSDIALOG_H

#include <QDialog>

namespace Ui 
{
    class EquilibriumAnalysisParametersDialog;
}

class EquilibriumAnalysisParametersDialog : public QDialog
{
    Q_OBJECT

public:

    // 
    // Constructor of the class.
    // 
    explicit EquilibriumAnalysisParametersDialog(QWidget *parent = nullptr);

    // 
    // Destructor of the class.
    // 
    ~EquilibriumAnalysisParametersDialog();

    // 
    // Returns the density value from the dialog.
    // @return double The density value from the dialog.
    // 
    double GetDensity() const;

    // 
    // Returns the friction coefficient value from the dialog.
    // @return double The friction coefficient from the dialog.
    // 
    double GetFrictionCoefficient() const;

    // 
    // Returns the length unit value from the dialog.
    // @return QString the length unit value from the dialog.
    // 
    QString GetLengthUnit() const;

    // 
    // Returns the compression weight value from the dialog.
    // @return double The compression weight value from the dialog.
    // 
    double GetCompressionWeight() const;

    // 
    // Returns the tension weight value from the dialog.
    // @return double The tension weight value from the dialog.
    // 
    double GetTensionWeight() const;

    // 
    // Returns the U tangential weight value from the dialog.
    // @return double The U tangential weight value from the dialog.
    // 
    double GetUTangentialWeight() const;

    // 
    // Returns the V tangential weight value from the dialog.
    // @return double The V tangential weight value from the dialog.
    // 
    double GetVTangentialWeight() const;

    // 
    // @return bool
    // 
    bool GetNormalize() const;

    // 
    // @return bool
    // 
    bool GetVerbose() const;

    // 
    // @return bool
    // 
    bool GetFiles() const;

private:

    Ui::EquilibriumAnalysisParametersDialog *ui;
};

#endif // EQUILIBRIUMANALYSISPARAMETERSDIALOG_H
