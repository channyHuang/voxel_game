#include <QApplication>

#include <QSurfaceFormat>

#include "mainwidget.h"

int main(int argc, char *argv[]) {
    QApplication a(argc, argv);

    QSurfaceFormat glFormat;
    glFormat.setVersion(3, 3);
    glFormat.setProfile(QSurfaceFormat::CoreProfile);
    QSurfaceFormat::setDefaultFormat(glFormat);

    MainWidget w;
    w.show();
    return a.exec();
}
