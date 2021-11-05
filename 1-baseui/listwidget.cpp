#include "listwidget.h"

ListWidget::ListWidget(QWidget *parent) : QListWidget(parent)
{

}

void ListWidget::addItem(QString qsText) {
    QListWidgetItem *item = new QListWidgetItem(this);
    item->setSizeHint(QSize(300, 50));

    ListItem *cusItem = new ListItem(this, ":/resources/cloud", qsText, ":/resources/upload", ":/resources/folder_small", ":/resources/world_small");
    this->setItemWidget(item, cusItem);
}
