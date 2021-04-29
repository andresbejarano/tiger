#ifndef TRUNCATEDCONEDIALOG_H
#define TRUNCATEDCONEDIALOG_H

#include <QDialog>

namespace Ui {
class TruncatedConeParametersDialog;
}

class TruncatedConeParametersDialog : public QDialog
{
    Q_OBJECT

public:

    /*
    Constructor of the class.
    */
    explicit TruncatedConeParametersDialog(QWidget *parent = nullptr);

    /*
    Destructor of the class.
    */
    ~TruncatedConeParametersDialog();

    /*
    Returns the selected axis value.
    @return QString The selected axis value.
    */
    QString GetAxis() const;

    /*
    Returns the bottom radius of the truncated cone.
    @return double The bottom radius of the truncated cone.
    */
    double GetBottomRadius() const;

    /*
    Returns the length of the truncated cone.
    @return double The length of the truncated cone.
    */
    double GetLength() const;

    /*
    Returns the number of length segments of the truncated cone.
    @return size_t The number of length segments of the truncated cone.
    */
    size_t GetLengthSegments() const;

    /*
    Returns the number of radial segments of the truncated cone.
    @return size_t The number of radial segments of the truncated cone.
    */
    size_t GetRadialSegments() const;

    /*
    Returns the top radius of the truncated cone.
    @return double The top radius of the truncated cone.
    */
    double GetTopRadius() const;

private:

    Ui::TruncatedConeParametersDialog *ui;
};

#endif // TRUNCATEDCONEDIALOG_H
