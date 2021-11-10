#ifndef LISTVIEW_H
#define LISTVIEW_H

#include <QWidget>
#include <QListView>

class ListView : public QListView
{
    Q_OBJECT
public:
    explicit ListView(QWidget *parent = nullptr);

signals:

};

#endif // LISTVIEW_H
