#ifndef MAINWIDGET_H
#define MAINWIDGET_H

#include <QWidget>
#include <QLabel>
#include <QTimer>
#include <QThread>

#include "renderControl/glwidget.h"
#include "terraingenerator_roblox.h"

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
    TerrainGenerator_Roblox::BiomesParam m_stBiomeParams;
};

#endif // MAINWIDGET_H
