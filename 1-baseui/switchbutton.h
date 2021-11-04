#ifndef SWITCHBUTTON_H
#define SWITCHBUTTON_H

#include <QWidget>
#include <QHBoxLayout>
#include <QVBoxLayout>
#include <QPainter>
#include <QLabel>
#include <QPalette>
#include <QPushButton>

class SwitchIcon : public QWidget {
    Q_OBJECT
public:
    SwitchIcon(QWidget *parent = nullptr, QColor background_color = Qt::gray, QColor forground_color = Qt::white, QColor checked_color = Qt::green);
    ~SwitchIcon() {}

    bool getChecked() { return bstate_; }

    void setBackground(QColor color);

signals:
    void beClicked();

protected:
    void paintEvent(QPaintEvent *event) override;
    void mousePressEvent(QMouseEvent *event) override;

private:
    QColor background_color_, forground_color_, checked_color_;
    bool bstate_ = false;
};

class SwitchButton : public QWidget {
    Q_OBJECT
public:
    enum LayoutType {
        BUTTON_LEFT,
        BUTTON_RIGHT,
        BUTTON_UP,
        BUTTON_DOWN
    };

    SwitchButton(QWidget *parent, QString qsLabel, LayoutType eLayoutType, QColor background_color = Qt::gray, QColor forground_color = Qt::white, QColor checked_color = Qt::green);
    ~SwitchButton();

    void setBackground(QColor color);
    void setTextColor(QColor color);
    void setSize(QSize size);

private:
    QBoxLayout *mainLayout = nullptr;

    QLabel *label = nullptr;
    SwitchIcon *button = nullptr;

    bool bstate = false;

};

#endif // SWITCHBUTTON_H
