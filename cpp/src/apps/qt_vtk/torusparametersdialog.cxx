#include "torusparametersdialog.h"
#include "ui_torusparametersdialog.h"

TorusParametersDialog::TorusParametersDialog(QWidget *parent) : QDialog(parent)
{
    setupUi();
}

TorusParametersDialog::~TorusParametersDialog()
{
    //delete ui;
}

QString TorusParametersDialog::getAxis() const
{
    return axisComboBox->currentText();
}

double TorusParametersDialog::getMajorRadius() const
{
    return majorRadiusDoubleSpinBox->value();
}

size_t TorusParametersDialog::getMajorRadialSegments() const
{
    return majorRadialSegmentsSpinBox->value();
}

double TorusParametersDialog::getMinorRadius() const
{
    return minorRadiusDoubleSpinBox->value();
}

size_t TorusParametersDialog::getMinorRadialSegments() const
{
    return minorRadialSegmentsSpinBox->value();
}

void TorusParametersDialog::setupUi() 
{
    // 
    this->setObjectName(QString::fromUtf8("TorusParametersDialog"));
    this->setWindowModality(Qt::ApplicationModal);
    this->resize(400, 208);
    this->setMinimumSize(QSize(400, 208));
    this->setWindowTitle(QString::fromUtf8("Torus Parameters"));

    // 
    gridLayout = new QGridLayout(this);
    gridLayout->setObjectName(QString::fromUtf8("gridLayout"));

    // 
    formLayout = new QFormLayout();
    formLayout->setObjectName(QString::fromUtf8("formLayout"));

    // Define the fields for the parameters
    majorRadiusLabel = new QLabel(this);
    majorRadiusLabel->setObjectName(QString::fromUtf8("majorRadiusLabel"));
    majorRadiusLabel->setText(QString::fromUtf8("Major Radius"));
    majorRadiusDoubleSpinBox = new QDoubleSpinBox(this);
    majorRadiusDoubleSpinBox->setObjectName(QString::fromUtf8("majorRadiusDoubleSpinBox"));
    majorRadiusDoubleSpinBox->setSingleStep(0.01);
    majorRadiusDoubleSpinBox->setValue(3.0);
    formLayout->setWidget(0, QFormLayout::LabelRole, majorRadiusLabel);
    formLayout->setWidget(0, QFormLayout::FieldRole, majorRadiusDoubleSpinBox);
    
    minorRadiusLabel = new QLabel(this);
    minorRadiusLabel->setObjectName(QString::fromUtf8("minorRadiusLabel"));
    minorRadiusLabel->setText(QString::fromUtf8("Minor Radius"));
    minorRadiusDoubleSpinBox = new QDoubleSpinBox(this);
    minorRadiusDoubleSpinBox->setObjectName(QString::fromUtf8("minorRadiusDoubleSpinBox"));
    minorRadiusDoubleSpinBox->setSingleStep(0.01);
    minorRadiusDoubleSpinBox->setValue(1.0);
    formLayout->setWidget(1, QFormLayout::LabelRole, minorRadiusLabel);
    formLayout->setWidget(1, QFormLayout::FieldRole, minorRadiusDoubleSpinBox);

    majorRadialSegmentsLabel = new QLabel(this);
    majorRadialSegmentsLabel->setObjectName(QString::fromUtf8("majorRadialSegmentsLabel"));
    majorRadialSegmentsLabel->setText(QString::fromUtf8("Major Radial Segments"));
    majorRadialSegmentsSpinBox = new QSpinBox(this);
    majorRadialSegmentsSpinBox->setObjectName(QString::fromUtf8("majorRadialSegmentsSpinBox"));
    majorRadialSegmentsSpinBox->setMinimum(1);
    majorRadialSegmentsSpinBox->setValue(3);
    formLayout->setWidget(2, QFormLayout::LabelRole, majorRadialSegmentsLabel);
    formLayout->setWidget(2, QFormLayout::FieldRole, majorRadialSegmentsSpinBox);
    
    minorRadialSegmentsLabel = new QLabel(this);
    minorRadialSegmentsLabel->setObjectName(QString::fromUtf8("minorRadialSegmentsLabel"));
    minorRadialSegmentsLabel->setText(QString::fromUtf8("Minor Radial Segments"));
    minorRadialSegmentsSpinBox = new QSpinBox(this);
    minorRadialSegmentsSpinBox->setObjectName(QString::fromUtf8("minorRadialSegmentsSpinBox"));
    minorRadialSegmentsSpinBox->setMinimum(1);
    minorRadialSegmentsSpinBox->setValue(3);
    formLayout->setWidget(3, QFormLayout::LabelRole, minorRadialSegmentsLabel);
    formLayout->setWidget(3, QFormLayout::FieldRole, minorRadialSegmentsSpinBox);
    
    axisLabel = new QLabel(this);
    axisLabel->setObjectName(QString::fromUtf8("axisLabel"));
    axisLabel->setText(QString::fromUtf8("Axis"));
    axisComboBox = new QComboBox(this);
    axisComboBox->addItem(QString::fromUtf8("X"));
    axisComboBox->addItem(QString::fromUtf8("Y"));
    axisComboBox->addItem(QString::fromUtf8("Z"));
    axisComboBox->setObjectName(QString::fromUtf8("axisComboBox"));
    formLayout->setWidget(4, QFormLayout::LabelRole, axisLabel);
    formLayout->setWidget(4, QFormLayout::FieldRole, axisComboBox);
    
    gridLayout->addLayout(formLayout, 0, 0, 1, 1);
    
    buttonBox = new QDialogButtonBox(this);
    buttonBox->setObjectName(QString::fromUtf8("buttonBox"));
    buttonBox->setOrientation(Qt::Horizontal);
    buttonBox->setStandardButtons(QDialogButtonBox::Cancel | QDialogButtonBox::Ok);
    
    gridLayout->addWidget(buttonBox, 1, 0, 1, 1);
    
    QWidget::setTabOrder(majorRadiusDoubleSpinBox, minorRadiusDoubleSpinBox);
    QWidget::setTabOrder(minorRadiusDoubleSpinBox, majorRadialSegmentsSpinBox);
    QWidget::setTabOrder(majorRadialSegmentsSpinBox, minorRadialSegmentsSpinBox);
    QWidget::setTabOrder(minorRadialSegmentsSpinBox, axisComboBox);
    
    QObject::connect(buttonBox, SIGNAL(accepted()), this, SLOT(accept()));
    QObject::connect(buttonBox, SIGNAL(rejected()), this, SLOT(reject()));
    
    QMetaObject::connectSlotsByName(this);
    this->setMinimumSize(this->size());
}
