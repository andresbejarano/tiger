#ifndef TILEOFFSETCLIPPINGPARAMETERSDIALOG_H
#define TILEOFFSETCLIPPINGPARAMETERSDIALOG_H

#include <QDialog>

namespace Ui {
class TileOffsetClippingParametersDialog;
}

class TileOffsetClippingParametersDialog : public QDialog
{
    Q_OBJECT

public:

    explicit TileOffsetClippingParametersDialog(QWidget *parent = nullptr);

    ~TileOffsetClippingParametersDialog();

    double GetExtrados() const;

    double GetIntrados() const;

private:
    Ui::TileOffsetClippingParametersDialog *ui;
};

#endif // TILEOFFSETCLIPPINGPARAMETERSDIALOG_H
