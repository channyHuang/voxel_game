#include "widget.h"

#include "popupdialog.h"
#include "switchbutton.h"

#include <QHBoxLayout>
#include <QVBoxLayout>
#include <QDebug>

Widget::Widget(QWidget *parent)
    : QWidget(parent)
{
    initWidget();

}

Widget::~Widget()
{
}

void Widget::testConfirmButton() {
    QString qsText = "xxxxxxxxxxxxxxxxxxxxxxxx";
    QVBoxLayout *mainLayout = new QVBoxLayout(this);

    ConfirmButton *button1 = new ConfirmButton(this, QString(":/resources/cloud"), qsText, ConfirmButton::TextUp);
    button1->setRoundCorner(15);
    button1->setBackground(Qt::gray);
    button1->setTextColor(Qt::red);
    QFont font1("Helvetica", 20);
    button1->setFont(font1);

    connect(button1, &ConfirmButton::bePress, [=]() {
        qDebug() << "be clicked";
    });

    ConfirmButton *button2 = new ConfirmButton(this, QString(":/resources/cloud"), qsText, ConfirmButton::TextDown);
    button2->setRoundCorner(15);
    button2->setBackground(Qt::gray);
    button2->setTextColor(Qt::red);
    QFont font2("Helvetica", 20);
    button2->setFont(font2);

    ConfirmButton *button3 = new ConfirmButton(this, QString(":/resources/cloud"), qsText, ConfirmButton::TextLeft);
    button3->setRoundCorner(15);
    button3->setBackground(Qt::gray);
    button3->setTextColor(Qt::red);
    QFont font3("Helvetica", 20);
    button3->setFont(font3);

    ConfirmButton *button4 = new ConfirmButton(this, QString(":/resources/cloud"), qsText, ConfirmButton::TextRight);
    button4->setRoundCorner(15);
    button4->setBackground(Qt::gray);
    button4->setTextColor(Qt::red);
    QFont font4("Helvetica", 20);
    button4->setFont(font4);

    mainLayout->addWidget(button1);
    mainLayout->addWidget(button2);
    mainLayout->addWidget(button3);
    mainLayout->addWidget(button4);

    setLayout(mainLayout);
}

void Widget::testPopupDialog() {
    QString qsMainText = "main text", qsSubText = "subtext";
    QString qsPosImage = ":/resources/upload", qsNegImage = ":/resources/cloud";

    PopupDialog *dialog = new PopupDialog(this, qsMainText, qsSubText, qsPosImage, qsNegImage);
    QFont font1("Helvetica", 25);
    dialog->setMainTextFont(font1);
    dialog->setMainTextColor(Qt::red);
    QFont font2("Helvetica", 23);
    dialog->setSubTextFont(font2);
    dialog->setSubTextColor(Qt::green);

    connect(dialog, &PopupDialog::beClickedPositive, [&](){
        qDebug() << "pos clicked";
    });
    connect(dialog, &PopupDialog::beClickedNegative, [&](){
        qDebug() << "neg clicked";
    });

    ConfirmButton *posButton = dialog->getConfirmButton();
    if (posButton != nullptr) {
        posButton->setRoundCorner(15);
        posButton->setBackground(Qt::gray);
        posButton->setTextColor(Qt::red);
        QFont font("Helvetica", 20);
        posButton->setFont(font);
    }

    QVBoxLayout *mainLayout = new QVBoxLayout(this);

    mainLayout->addWidget(dialog);
    setLayout(mainLayout);
}

void Widget::testSwitchButton() {
    SwitchButton *button = new SwitchButton(nullptr, "xxxxx", SwitchButton::BUTTON_RIGHT);

    QVBoxLayout *mainLayout = new QVBoxLayout(this);

    mainLayout->addWidget(button);

    setLayout(mainLayout);
}

void Widget::initWidget() {
    testSwitchButton();
}
