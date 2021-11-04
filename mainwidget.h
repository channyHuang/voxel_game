#ifndef MAINWIDGET_H
#define MAINWIDGET_H

#include <QWidget>
#include <QLabel>

#include "glwidget.h"

class MainWidget : public QWidget {
    Q_OBJECT
public:
    MainWidget(QWidget *parent = nullptr);
    ~MainWidget();

private:
    GlWidget *m_pGlWidget = nullptr;
};

#endif // MAINWIDGET_H
