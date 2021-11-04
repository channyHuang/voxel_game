#include "switchbutton.h"

SwitchIcon::SwitchIcon(QWidget *parent /* = nullptr */,
                            QColor background_color /* = Qt::gray */,
                            QColor forground_color /* = Qt::white */,  QColor checked_color /* = Qt::green */) : QWidget(parent) {
    Q_UNUSED(parent);
    background_color_ = background_color;
    forground_color_ = forground_color;
    checked_color_ = checked_color;
}

void SwitchIcon::paintEvent(QPaintEvent *event) {
    QPainter painter(this);
    painter.setRenderHint(QPainter::Antialiasing);
    painter.setBrush(QBrush(bstate_ ? checked_color_ : background_color_));
    painter.setPen(Qt::transparent);

    QRect rect = this->rect();
    if (rect.width() < rect.height()) {
        rect.setWidth(rect.height());
    }
    painter.drawRoundedRect(rect, rect.height() / 2, rect.height() / 2);

    int stPos = (bstate_ ? rect.width() - rect.height() : 0);
    QRect circle = QRect(stPos, 0, rect.height(), this->rect().height());
    painter.setBrush(QBrush(forground_color_));
    painter.drawRoundedRect(circle, rect.height() / 2, rect.height() / 2);

    QWidget::paintEvent(event);
}

void SwitchIcon::mousePressEvent(QMouseEvent *event) {
    bstate_ = !bstate_;
    QWidget::mousePressEvent(event);
    update();
}

void SwitchIcon::setBackground(QColor color) {
    QPalette palette = this->palette();
    palette.setColor(QPalette::Background, color);
    this->setPalette(palette);
}

SwitchButton::SwitchButton(QWidget *parent, QString qsLabel, LayoutType eLayoutType,
    QColor background_color, QColor forground_color, QColor checked_color)
    : QWidget(parent) {
    button = new SwitchIcon(this, background_color, forground_color, checked_color);
    button->setMinimumSize(QSize(button->height(), button->height()));
    if (qsLabel.length() > 0) {
        label = new QLabel(qsLabel);
    }

    switch (eLayoutType)
    {
    case SwitchButton::BUTTON_LEFT:
        mainLayout = new QHBoxLayout();
        mainLayout->addWidget(button);
        mainLayout->addWidget(label);
        break;
    case SwitchButton::BUTTON_RIGHT:
        mainLayout = new QHBoxLayout();
        mainLayout->addWidget(label);
        mainLayout->addWidget(button);
        break;
    case SwitchButton::BUTTON_UP:
        mainLayout = new QVBoxLayout();
        mainLayout->addWidget(button);
        mainLayout->addWidget(label);
        break;
    case SwitchButton::BUTTON_DOWN:
        mainLayout = new QVBoxLayout();
        mainLayout->addWidget(label);
        mainLayout->addWidget(button);
        break;
    default:
        break;
    }


    setLayout(mainLayout);
}

SwitchButton::~SwitchButton() {}

void SwitchButton::setTextColor(QColor color) {
    if (label == nullptr) return;
    QPalette palette = label->palette();
    palette.setColor(QPalette::Foreground, color);
    label->setPalette(palette);
}

void SwitchButton::setSize(QSize size) {
    this->setFixedSize(size);
}

void SwitchButton::setBackground(QColor color) {
    if (button != nullptr) {
        button->setBackground(color);
    }
}
