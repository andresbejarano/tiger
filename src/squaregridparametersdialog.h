#ifndef SQUAREGRIDPARAMETERSDIALOG_H
#define SQUAREGRIDPARAMETERSDIALOG_H

#include <QDialog>

namespace Ui {
class SquareGridParametersDialog;
}

class SquareGridParametersDialog : public QDialog
{
    Q_OBJECT

public:

    explicit SquareGridParametersDialog(QWidget *parent = nullptr);

    ~SquareGridParametersDialog();

    double GetWidth() const;

    double GetHeight() const;

    size_t GetWidthSegments() const;

    size_t GetHeightSegments() const;

    QString GetPlane() const;

private:

    Ui::SquareGridParametersDialog *ui;
};

#endif // SQUAREGRIDPARAMETERSDIALOG_H
