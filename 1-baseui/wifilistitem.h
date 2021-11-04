#ifndef WIFILISTITEM_H
#define WIFILISTITEM_H

#include <QWidget>
#include <QListWidgetItem>

class wifilistitem : public QListWidgetItem
{
    Q_OBJECT
public:
    explicit wifilistitem(QListWidget *parent = nullptr);

protected:

signals:

};

#endif // WIFILISTITEM_H
