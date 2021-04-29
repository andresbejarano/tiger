#ifndef TORUSPARAMETERSDIALOG_H
#define TORUSPARAMETERSDIALOG_H

#include <QComboBox>
#include <QDialog>
#include <QDialogButtonBox>
#include <QDoubleSpinBox>
#include <QFormLayout>
#include <QLabel>
#include <QSpinBox>

class TorusParametersDialog : public QDialog
{
    Q_OBJECT

private:

    QGridLayout* gridLayout;
    QFormLayout* formLayout;

    QLabel* majorRadiusLabel;
    QDoubleSpinBox* majorRadiusDoubleSpinBox;

    QLabel* majorRadialSegmentsLabel;
    QSpinBox* majorRadialSegmentsSpinBox;

    QLabel* minorRadiusLabel;
    QDoubleSpinBox* minorRadiusDoubleSpinBox;

    QLabel* minorRadialSegmentsLabel;
    QSpinBox* minorRadialSegmentsSpinBox;

    QLabel* axisLabel;
    QComboBox* axisComboBox;
    
    QDialogButtonBox* buttonBox;

public:

    /*
    Constructor of the class.
    */
    explicit TorusParametersDialog(QWidget *parent = nullptr);

    /*
    Destructor of the class.
    */
    ~TorusParametersDialog();

    /*
    Returns the selected axis value.
    @return QString The selected axis value.
    */
    QString getAxis() const;

    /*
    Returns the major radius of the torus.
    @return double The major radius of the torus.
    */
    double getMajorRadius() const;

    /// <summary>
    /// Returns the number of major radial segments of the torus.
    /// </summary>
    /// <returns>size_t The number of major radial segments of the torus.</returns>
    size_t getMajorRadialSegments() const;

    /// <summary>
    /// Returns the minor radius of the torus.
    /// </summary>
    /// <returns>double The minor radius of the torus.</returns>
    double getMinorRadius() const;

    /// <summary>
    /// Returns the number of minor radial segments of the torus.
    /// </summary>
    /// <returns>size_t The number of minor radial segments of the torus.</returns>
    size_t getMinorRadialSegments() const;

private:

    /// <summary>
    /// 
    /// </summary>
    void setupUi();
};

#endif // TORUSPARAMETERSDIALOG_H
