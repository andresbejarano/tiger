#ifndef SADDLEPARAMETERSDIALOG_H
#define SADDLEPARAMETERSDIALOG_H

#include <QDialog>

namespace Ui {
class SaddleParametersDialog;
}

class SaddleParametersDialog : public QDialog
{
    Q_OBJECT

public:

    explicit SaddleParametersDialog(QWidget *parent = nullptr);

    ~SaddleParametersDialog();

    double GetWidth() const;

    double GetHeight() const;

    size_t GetWidthSegments() const;

    size_t GetHeightSegments() const;

    QString GetPlane() const;

private:
    Ui::SaddleParametersDialog *ui;
};

#endif // SADDLEPARAMETERSDIALOG_H
