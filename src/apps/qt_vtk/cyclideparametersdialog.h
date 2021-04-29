#ifndef _CYCLIDEPARAMETERSDIALOG_H_
#define _CYCLIDEPARAMETERSDIALOG_H_

#include <QDialog>
#include <QDialogButtonBox>
#include <QDoubleSpinBox>
#include <QFormLayout>
#include <QLabel>
#include <QSpinBox>

class CyclideParametersDialog : public QDialog
{
    Q_OBJECT

private:

    QGridLayout* gridLayout;
    QFormLayout* formLayout;

    QLabel* aLabel;
    QDoubleSpinBox* aDoubleSpinBox;

    QLabel* bLabel;
    QDoubleSpinBox* bDoubleSpinBox;

    QLabel* cLabel;
    QDoubleSpinBox* cDoubleSpinBox;

    QLabel* dLabel;
    QDoubleSpinBox* dDoubleSpinBox;

    QLabel* majorRadialSegmentsLabel;
    QSpinBox* majorRadialSegmentsSpinBox;

    QLabel* minorRadialSegmentsLabel;
    QSpinBox* minorRadialSegmentsSpinBox;

    QDialogButtonBox* buttonBox;

public:

    /// <summary>
    /// Constructor of the class.
    /// </summary>
    /// <param name="parent"></param>
    explicit CyclideParametersDialog(QWidget *parent = nullptr);

    ///
    ~CyclideParametersDialog();

    /// <summary>
    /// 
    /// </summary>
    /// <returns></returns>
    double getA() const;

    /// <summary>
    /// 
    /// </summary>
    /// <returns></returns>
    double getB() const;

    /// <summary>
    /// 
    /// </summary>
    /// <returns></returns>
    double getC() const;

    /// <summary>
    /// 
    /// </summary>
    /// <returns></returns>
    double getD() const;

    /// <summary>
    /// 
    /// </summary>
    /// <returns></returns>
    size_t getMajorRadialSegments() const;

    /// <summary>
    /// 
    /// </summary>
    /// <returns></returns>
    size_t getMinorRadialSegments() const;

private:

    /// <summary>
    /// 
    /// </summary>
    void setupUi();
};

#endif // _CYCLIDEPARAMETERSDIALOG_H_
