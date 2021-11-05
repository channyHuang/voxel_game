#ifndef WIFILISTITEM_H
#define WIFILISTITEM_H

#include <QWidget>
#include <QHBoxLayout>
#include <QVBoxLayout>
#include <QLabel>
#include <QListWidget>

class ListItem : public QWidget
{
    Q_OBJECT
public:
    explicit ListItem(QListWidget *parent = nullptr, QString qsIcon = "", QString qsText = "", QString qsState = "", QString qsSignal = "", QString qsInfoIcon = "");
    ListItem(QWidget *parent, QStringList qslInputs);
protected:

signals:

};

#endif // WIFILISTITEM_H
