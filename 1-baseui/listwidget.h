#ifndef LISTWIDGET_H
#define LISTWIDGET_H

#include <QObject>
#include <QWidget>
#include <QListWidget>
#include <QListWidgetItem>

#include "listitem.h"

class ListWidget : public QListWidget
{
    Q_OBJECT
public:
    explicit ListWidget(QWidget *parent = nullptr);

    void addItem(QString qsText);
signals:

private:

};

#endif // LISTWIDGET_H
