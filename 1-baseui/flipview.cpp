#include "flipview.h"

FlipView::FlipView(QWidget *parent) : QWidget(parent)
{
    QLabel *textLabel = new QLabel("text in front of background image");
    QLabel *subTextLabel = new QLabel("sub text in front of image");

    background_label = new QLabel(this);
    background_label->setPixmap(QPixmap::fromImage(QImage(":/resources/500x500")));

    //QLayout *labelLayout = background_label->layout();
    //if (labelLayout == nullptr) {
    //    labelLayout = new QVBoxLayout;
    //}
    QHBoxLayout *labelLayout = new QHBoxLayout;

    QVBoxLayout *textLayout = new QVBoxLayout;
    textLayout->addStretch();
    textLayout->addWidget(textLabel);
    textLayout->addStretch();
    textLayout->addWidget(subTextLabel);
    textLayout->addStretch();

    labelLayout->addStretch();
    labelLayout->addLayout(textLayout);
    labelLayout->addStretch();
    background_label->setLayout(labelLayout);

    QVBoxLayout *mainLayout = new QVBoxLayout;
    mainLayout->addWidget(background_label);
    setLayout(mainLayout);
}
