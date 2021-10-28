#ifndef WIDGET_H
#define WIDGET_H

#include <QWidget>

class Widget : public QWidget
{
    Q_OBJECT

public:
    Widget(QWidget *parent = nullptr);
    ~Widget();

    void testConfirmButton();
    void testPopupDialog();
private:
    void initWidget();

private:
};
#endif // WIDGET_H
