#include "listitem.h"

ListItem::ListItem(QListWidget *parent, QString qsIcon, QString qsText, QString qsState, QString qsSignal, QString qsInfoIcon)
    : QWidget(parent)
{
    setMinimumHeight(50);

    QHBoxLayout *mainLayout = new QHBoxLayout;

    QLabel *iconLabel = new QLabel(this);
    iconLabel->setPixmap(QPixmap::fromImage(QImage(qsIcon)));

    QLabel *stateLabel = new QLabel(this);
    stateLabel->setPixmap(QPixmap::fromImage(QImage(qsState)));

    QLabel *signalLabel = new QLabel(this);
    signalLabel->setPixmap(QPixmap::fromImage(QImage(qsSignal)));

    QLabel *infoLabel = new QLabel(this);
    infoLabel->setPixmap(QPixmap::fromImage(QImage(qsInfoIcon)));


    mainLayout->addWidget(iconLabel);
    mainLayout->addWidget(new QLabel(qsText));
    mainLayout->addStretch();
    mainLayout->addWidget(stateLabel);
    mainLayout->addWidget(signalLabel);
    mainLayout->addWidget(infoLabel);

    setLayout(mainLayout);
}

ListItem::ListItem(QWidget *parent, QStringList qslInputs) {

}
