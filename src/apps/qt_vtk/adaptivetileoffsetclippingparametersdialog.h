#ifndef ADAPTIVETILEOFFSETCLIPPINGPARAMETERSDIALOG_H
#define ADAPTIVETILEOFFSETCLIPPINGPARAMETERSDIALOG_H

#include <QDialog>

namespace Ui {
class AdaptiveTileOffsetClippingParametersDialog;
}

class AdaptiveTileOffsetClippingParametersDialog : public QDialog
{
    Q_OBJECT

public:

    /*
    */
    explicit AdaptiveTileOffsetClippingParametersDialog(QWidget *parent = nullptr);

    /*
    */
    ~AdaptiveTileOffsetClippingParametersDialog();

    /*
    */
    QString GetBottomFunction() const;

    /*
    */
    QString GetTopFunction() const;

private:
    Ui::AdaptiveTileOffsetClippingParametersDialog *ui;
};

#endif // ADAPTIVETILEOFFSETCLIPPINGPARAMETERSDIALOG_H
