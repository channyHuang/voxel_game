#ifndef CONFIRMBUTTON_H
#define CONFIRMBUTTON_H

#include <QVBoxLayout>
#include <QHBoxLayout>
#include <QLabel>
#include <QPixmap>
#include <QString>
#include <QPalette>
#include <QPainter>
#include <QSizePolicy>

class ConfirmButton : public QWidget
{
    Q_OBJECT
public:
    ConfirmButton(QWidget *parent = nullptr);

    enum TextDirection {
        TextUp,
        TextDown,
        TextLeft,
        TextRight
    };

    ConfirmButton(QWidget *parent = nullptr, QString qsImageName = "", QString qsText = "", TextDirection eTextDirection = TextRight);
    // background
    void setBackground(QColor color);
    void setRoundCorner(qreal radius);
    // text
    void setTextColor(QColor color);
    void setTextFont(QFont font);
    // icon
    void setIconSize(QSize size);

signals:
    void bePress();

protected:
    virtual void mousePressEvent(QMouseEvent *event) override;
    virtual void mouseMoveEvent(QMouseEvent *event) override;
    virtual void mouseReleaseEvent(QMouseEvent *event) override;
    virtual void paintEvent(QPaintEvent *event) override;

private:
    void initWidget();

private:
    QLabel *m_textLabel = nullptr;
    QLabel *m_iconLabel = nullptr;
    TextDirection m_eTextDirection = TextDown;

    QBoxLayout *mainLayout = nullptr;
    QBoxLayout *textLayout = nullptr, *iconLayout = nullptr;
    qreal m_xRadius = 0, m_yRadius = 0;
};

#endif // CONFIRMBUTTON_H
