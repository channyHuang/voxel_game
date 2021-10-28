#include "confirmbutton.h"

ConfirmButton::ConfirmButton(QWidget *parent)
    : QWidget(parent) {

}

ConfirmButton::ConfirmButton(QWidget *parent, QString qsImageName, QString qsText, TextDirection eTextDirection)
    : QWidget(parent) {
    setAttribute(Qt::WA_TranslucentBackground);

    m_eTextDirection = eTextDirection;
    if (qsImageName.length() > 0) {
        m_iconLabel = new QLabel();
        QImage iconImage = QImage(qsImageName);
         m_iconLabel->setPixmap(QPixmap::fromImage(iconImage));
         m_iconLabel->setSizePolicy(QSizePolicy::Preferred, QSizePolicy::Preferred);
    }
    if (qsText.length() > 0) {
        m_textLabel = new QLabel(qsText);
        m_textLabel->setSizePolicy(QSizePolicy::Preferred, QSizePolicy::Preferred);
    }

    initWidget();
}

void ConfirmButton::initWidget() {
    if (m_eTextDirection == TextUp || m_eTextDirection == TextDown) {
        mainLayout = new QVBoxLayout(this);
    } else {
        mainLayout = new QHBoxLayout(this);
    }

    switch(m_eTextDirection) {
    case TextUp:
        textLayout = new QHBoxLayout();
        iconLayout = new QHBoxLayout();
        if (m_textLabel != nullptr) {
            textLayout->addStretch();
            textLayout->addWidget(m_textLabel);
            textLayout->addStretch();
        }
        mainLayout->addLayout(textLayout);
        mainLayout->addStretch();
        if (m_iconLabel != nullptr) {
            iconLayout->addStretch();
            iconLayout->addWidget(m_iconLabel);
            iconLayout->addStretch();
        }
        mainLayout->addLayout(iconLayout);
        break;
    case TextLeft:
        textLayout = new QVBoxLayout();
        iconLayout = new QVBoxLayout();
        if (m_textLabel != nullptr) {
            textLayout->addStretch();
            textLayout->addWidget(m_textLabel);
            textLayout->addStretch();
        }
        mainLayout->addLayout(textLayout, Qt::AlignCenter);
        mainLayout->addStretch();
        if (m_iconLabel != nullptr) {
            iconLayout->addStretch();
            iconLayout->addWidget(m_iconLabel, Qt::AlignCenter);
            iconLayout->addStretch();
        }
        mainLayout->addLayout(iconLayout);
        break;
    case TextDown:
        iconLayout = new QHBoxLayout();
        textLayout = new QHBoxLayout();
        if (m_iconLabel != nullptr) {
            iconLayout->addStretch();
            iconLayout->addWidget(m_iconLabel);
            iconLayout->addStretch();
        }
        mainLayout->addLayout(iconLayout);
        mainLayout->addStretch();
        if (m_textLabel != nullptr) {
            textLayout->addStretch();
            textLayout->addWidget(m_textLabel);
            textLayout->addStretch();
        }
        mainLayout->addLayout(textLayout);
        break;
    case TextRight:
        textLayout = new QVBoxLayout();
        iconLayout = new QVBoxLayout();
        if (m_iconLabel != nullptr) {
            iconLayout->addStretch();
            iconLayout->addWidget(m_iconLabel, Qt::AlignCenter);
            iconLayout->addStretch();
        }
        mainLayout->addLayout(iconLayout);
        mainLayout->addStretch();
        if (m_textLabel != nullptr) {
            textLayout->addStretch();
            textLayout->addWidget(m_textLabel);
            textLayout->addStretch();
        }
        mainLayout->addLayout(textLayout);
        break;
    default:
        break;
    }

    setLayout(mainLayout);
}

void ConfirmButton::setBackground(QColor color) {
    QPalette palette = this->palette();
    palette.setColor(QPalette::Background, color);
    setPalette(palette);
}

void ConfirmButton::setRoundCorner(qreal radius) {
    m_xRadius = m_yRadius = radius;
}

void ConfirmButton::setTextColor(QColor color) {
    QPalette palette = this->palette();
    palette.setColor(QPalette::Foreground, color);
    setPalette(palette);
}

void ConfirmButton::setTextFont(QFont font) {
    if (m_textLabel != nullptr) {
        m_textLabel->setFont(font);
    }
}

void ConfirmButton::setIconSize(QSize size) {
    if (m_iconLabel != nullptr) {
        m_iconLabel->setMaximumSize(size);
    }
}

void ConfirmButton::paintEvent(QPaintEvent *event) {
    QPalette palette = this->palette();

    QPainter painter(this);
    painter.setRenderHint(QPainter::Antialiasing);
    painter.setBrush(QBrush(palette.background().color()));
    painter.setPen(Qt::transparent);

    QRect rect = this->rect();
    painter.drawRoundedRect(rect, m_xRadius, m_yRadius);

    QWidget::paintEvent(event);
}

void ConfirmButton::mousePressEvent(QMouseEvent *event) {
    QWidget::mousePressEvent(event);
    emit bePress();
}

void ConfirmButton::mouseMoveEvent(QMouseEvent *event) {
    QWidget::mouseMoveEvent(event);
}

void ConfirmButton::mouseReleaseEvent(QMouseEvent *event) {
    QWidget::mouseReleaseEvent(event);
}
