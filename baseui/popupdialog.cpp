#include "popupdialog.h"

PopupDialog::PopupDialog(QWidget *parent) : QWidget(parent)
{

}

PopupDialog::PopupDialog(QWidget *parent,
                         QString qsMainText, QString qsSubText,
                         QString qsPositiveIcon, QString qsNegativeIcon,
                         QString qsPositiveText, QString qsNegativeText)
    : QWidget(parent) {

    if (qsMainText.length() > 0) {
        mainTextLabel = new QLabel(qsMainText);
    }
    if (qsSubText.length() > 0) {
        subTextLabel = new QLabel(qsSubText);
    }
    if (qsPositiveIcon.length() > 0 || qsPositiveText.length() > 0) {
        posButton = new ConfirmButton(this, qsPositiveIcon, qsPositiveText, ConfirmButton::TextDown);
        connect(posButton, &ConfirmButton::bePress, this, &PopupDialog::beClickedPositive);
    }
    if (qsNegativeIcon.length() > 0 || qsNegativeText.length() > 0) {
        negButton = new ConfirmButton(this, qsNegativeIcon, qsNegativeText, ConfirmButton::TextDown);
        connect(negButton, &ConfirmButton::bePress, this, &PopupDialog::beClickedNegative);
    }

    initWidget();
}

PopupDialog::~PopupDialog() {}

void PopupDialog::initWidget() {
    mainLayout = new QVBoxLayout(this);

    if (mainTextLabel != nullptr) {
        QHBoxLayout *mainTextLayout = new QHBoxLayout();
        mainTextLayout->addStretch();
        mainTextLayout->addWidget(mainTextLabel);
        mainTextLayout->addStretch();

        mainLayout->addLayout(mainTextLayout);
    }
    if (subTextLabel != nullptr) {
        QHBoxLayout *subTextLayout = new QHBoxLayout();
        subTextLayout->addStretch();
        subTextLayout->addWidget(subTextLabel);
        subTextLayout->addStretch();

        mainLayout->addLayout(subTextLayout);
    }

    QHBoxLayout *buttonLayout = new QHBoxLayout();
    buttonLayout->addWidget(posButton);
    buttonLayout->addWidget(negButton);
    mainLayout->addStretch();
    mainLayout->addLayout(buttonLayout);
    mainLayout->addStretch();

    setLayout(mainLayout);
}

void PopupDialog::paintEvent(QPaintEvent *event) {
    QPalette palette = this->palette();

    QPainter painter(this);
    painter.setRenderHint(QPainter::Antialiasing);
    painter.setBrush(QBrush(palette.background().color()));
    painter.setPen(Qt::transparent);

    QRect rect = this->rect();
    painter.drawRoundedRect(rect, m_xRadius, m_yRadius);

    QWidget::paintEvent(event);
}

void PopupDialog::setMainTextColor(QColor color) {
    if (mainTextLabel != nullptr) {
        QPalette palette = mainTextLabel->palette();
        palette.setColor(QPalette::Foreground, color);
        mainTextLabel->setPalette(palette);
    }
}

void PopupDialog::setMainTextFont(QFont font) {
    if (mainTextLabel != nullptr) {
        mainTextLabel->setFont(font);
    }
}

void PopupDialog::setSubTextColor(QColor color) {
    if (subTextLabel != nullptr) {
        QPalette palette = subTextLabel->palette();
        palette.setColor(QPalette::Foreground, color);
        subTextLabel->setPalette(palette);
    }
}

void PopupDialog::setSubTextFont(QFont font) {
    if (subTextLabel != nullptr) {
        subTextLabel->setFont(font);
    }
}

void PopupDialog::setBackground(QColor color) {
    QPalette palette = this->palette();
    palette.setColor(QPalette::Background, color);
    setPalette(palette);
}

void PopupDialog::setRoundCorner(qreal radius) {
    m_xRadius = m_yRadius = radius;
}

ConfirmButton* PopupDialog::getConfirmButton(bool bPositive) {
    return (bPositive ? posButton : negButton);
}
