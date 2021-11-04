#ifndef POPUPDIALOG_H
#define POPUPDIALOG_H

#include <QWidget>
#include <QHBoxLayout>
#include <QVBoxLayout>
#include <QLabel>

#include "confirmbutton.h"

class PopupDialog : public QWidget
{
    Q_OBJECT
public:
    explicit PopupDialog(QWidget *parent = nullptr);
    PopupDialog(QWidget *parent = nullptr,
                QString qsMainText = "", QString qsSubText = "",
                QString qsPositiveIcon = "", QString qsnegativeIcon = "",
                QString qsPositiveText = tr("Confirm"), QString qsNegativeText = tr("Cancel"));
    ~PopupDialog();

    // background
    void setBackground(QColor color);
    void setRoundCorner(qreal radius);
    // text
    void setMainTextColor(QColor color);
    void setMainTextFont(QFont font);
    void setSubTextColor(QColor color);
    void setSubTextFont(QFont font);

    ConfirmButton* getConfirmButton(bool bPositive = true);

signals:
    void beClickedPositive();
    void beClickedNegative();

protected:
    virtual void paintEvent(QPaintEvent *event) override;

private:
    void initWidget();

private:
    QVBoxLayout *mainLayout = nullptr;
    QLabel *posIconLabel = nullptr, *negIconLabel = nullptr;
    QLabel *mainTextLabel = nullptr, *subTextLabel = nullptr;
    ConfirmButton *posButton = nullptr, *negButton = nullptr;
    qreal m_xRadius = 0, m_yRadius = 0;
};

#endif // POPUPDIALOG_H
