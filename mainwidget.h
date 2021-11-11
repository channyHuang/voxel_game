#ifndef MAINWIDGET_H
#define MAINWIDGET_H

#include <QWidget>
#include <QLabel>
#include <QTimer>
#include <QThread>

#include "glwidget.h"

class MainWidget : public QWidget {
    Q_OBJECT
public:
    MainWidget(QWidget *parent = nullptr);
    ~MainWidget();

    bool init();
    void initWidgets();

signals:
    void notice(int event);

private:
    GlWidget *m_pGlWidget = nullptr;
    QTimer timer;
    QThread *terrain_thread = nullptr;
};

#endif // MAINWIDGET_H
