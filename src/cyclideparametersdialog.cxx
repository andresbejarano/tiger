#include "cyclideparametersdialog.h"

CyclideParametersDialog::CyclideParametersDialog(QWidget* parent) : QDialog(parent)
{
    setupUi();
}

CyclideParametersDialog::~CyclideParametersDialog() 
{
    delete aLabel;
    delete aDoubleSpinBox;
    delete bLabel;
    delete bDoubleSpinBox;
    delete cLabel;
    delete cDoubleSpinBox;
    delete dLabel;
    delete dDoubleSpinBox;
    delete majorRadialSegmentsLabel;
    delete majorRadialSegmentsSpinBox;
    delete minorRadialSegmentsLabel;
    delete minorRadialSegmentsSpinBox;
    delete buttonBox;
    delete formLayout;
    delete gridLayout;
}

void CyclideParametersDialog::setupUi() 
{
    // 
    this->setObjectName(QString::fromUtf8("CyclideParametersDialog"));
    this->setWindowModality(Qt::ApplicationModal);
    this->resize(400, 208);
    this->setMinimumSize(QSize(400, 208));
    this->setWindowTitle(QString::fromUtf8("Cyclide Parameters"));

    // 
    gridLayout = new QGridLayout(this);
    gridLayout->setObjectName(QString::fromUtf8("gridLayout"));

    // 
    formLayout = new QFormLayout();
    formLayout->setObjectName(QString::fromUtf8("formLayout"));

    // Define the fields for the a parameter
    aLabel = new QLabel(this);
    aLabel->setObjectName(QString::fromUtf8("aLabel"));
    aLabel->setText(QString::fromUtf8("a"));
    aDoubleSpinBox = new QDoubleSpinBox(this);
    aDoubleSpinBox->setObjectName(QString::fromUtf8("aDoubleSpinBox"));
    aDoubleSpinBox->setSingleStep(0.01);
    aDoubleSpinBox->setValue(1);
    formLayout->setWidget(0, QFormLayout::LabelRole, aLabel);
    formLayout->setWidget(0, QFormLayout::FieldRole, aDoubleSpinBox);

    bLabel = new QLabel(this);
    bLabel->setObjectName(QString::fromUtf8("bLabel"));
    bLabel->setText(QString::fromUtf8("b"));
    bDoubleSpinBox = new QDoubleSpinBox(this);
    bDoubleSpinBox->setObjectName(QString::fromUtf8("bDoubleSpinBox"));
    bDoubleSpinBox->setSingleStep(0.01);
    bDoubleSpinBox->setValue(1);
    formLayout->setWidget(1, QFormLayout::LabelRole, bLabel);
    formLayout->setWidget(1, QFormLayout::FieldRole, bDoubleSpinBox);

    cLabel = new QLabel(this);
    cLabel->setObjectName(QString::fromUtf8("cLabel"));
    cLabel->setText(QString::fromUtf8("c"));
    cDoubleSpinBox = new QDoubleSpinBox(this);
    cDoubleSpinBox->setObjectName(QString::fromUtf8("cDoubleSpinBox"));
    cDoubleSpinBox->setSingleStep(0.01);
    cDoubleSpinBox->setValue(0.2);
    formLayout->setWidget(2, QFormLayout::LabelRole, cLabel);
    formLayout->setWidget(2, QFormLayout::FieldRole, cDoubleSpinBox);

    dLabel = new QLabel(this);
    dLabel->setObjectName(QString::fromUtf8("dLabel"));
    dLabel->setText(QString::fromUtf8("d"));
    dDoubleSpinBox = new QDoubleSpinBox(this);
    dDoubleSpinBox->setObjectName(QString::fromUtf8("dDoubleSpinBox"));
    dDoubleSpinBox->setSingleStep(0.01);
    dDoubleSpinBox->setValue(0.4);
    formLayout->setWidget(3, QFormLayout::LabelRole, dLabel);
    formLayout->setWidget(3, QFormLayout::FieldRole, dDoubleSpinBox);

    majorRadialSegmentsLabel = new QLabel(this);
    majorRadialSegmentsLabel->setObjectName(QString::fromUtf8("majorRadialSegmentsLabel"));
    majorRadialSegmentsLabel->setText(QString::fromUtf8("Major radial segments"));
    majorRadialSegmentsSpinBox = new QSpinBox(this);
    majorRadialSegmentsSpinBox->setObjectName(QString::fromUtf8("majorRadialSegmentsSpinBox"));
    majorRadialSegmentsSpinBox->setValue(30);
    formLayout->setWidget(4, QFormLayout::LabelRole, majorRadialSegmentsLabel);
    formLayout->setWidget(4, QFormLayout::FieldRole, majorRadialSegmentsSpinBox);

    minorRadialSegmentsLabel = new QLabel(this);
    minorRadialSegmentsLabel->setObjectName(QString::fromUtf8("minorRadialSegmentsLabel"));
    minorRadialSegmentsLabel->setText(QString::fromUtf8("Minor radial segments"));
    minorRadialSegmentsSpinBox = new QSpinBox(this);
    minorRadialSegmentsSpinBox->setObjectName(QString::fromUtf8("minorRadialSegmentsSpinBox"));
    minorRadialSegmentsSpinBox->setValue(20);
    formLayout->setWidget(5, QFormLayout::LabelRole, minorRadialSegmentsLabel);
    formLayout->setWidget(5, QFormLayout::FieldRole, minorRadialSegmentsSpinBox);


    //...

    gridLayout->addLayout(formLayout, 0, 0, 1, 1);

    buttonBox = new QDialogButtonBox(this);
    buttonBox->setObjectName(QString::fromUtf8("buttonBox"));
    buttonBox->setOrientation(Qt::Horizontal);
    buttonBox->setStandardButtons(QDialogButtonBox::Cancel | QDialogButtonBox::Ok);

    gridLayout->addWidget(buttonBox, 1, 0, 1, 1);

    QWidget::setTabOrder(aDoubleSpinBox, bDoubleSpinBox);
    QWidget::setTabOrder(bDoubleSpinBox, cDoubleSpinBox);
    QWidget::setTabOrder(cDoubleSpinBox, dDoubleSpinBox);
    QWidget::setTabOrder(dDoubleSpinBox, majorRadialSegmentsSpinBox);
    QWidget::setTabOrder(majorRadialSegmentsSpinBox, minorRadialSegmentsSpinBox);

    QObject::connect(buttonBox, SIGNAL(accepted()), this, SLOT(accept()));
    QObject::connect(buttonBox, SIGNAL(rejected()), this, SLOT(reject()));

    QMetaObject::connectSlotsByName(this);

    this->setMinimumSize(this->size());
}

double CyclideParametersDialog::getA() const 
{
    return this->aDoubleSpinBox->value();
}

double CyclideParametersDialog::getB() const 
{
    return this->bDoubleSpinBox->value();
}

double CyclideParametersDialog::getC() const
{
    return cDoubleSpinBox->value();
}

double CyclideParametersDialog::getD() const 
{
    return dDoubleSpinBox->value();
}

size_t CyclideParametersDialog::getMajorRadialSegments() const 
{
    return majorRadialSegmentsSpinBox->value();
}

size_t CyclideParametersDialog::getMinorRadialSegments() const 
{
    return minorRadialSegmentsSpinBox->value();
}