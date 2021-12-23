#ifndef OBJWORKER_H
#define OBJWORKER_H

#include <vector>
#include <QDir>
#include <QString>
#include <QDebug>
#include <QFile>
#include <QFileInfo>
#include <fstream>

class ObjWorker
{
public:
    ObjWorker();
    ~ObjWorker();

    static bool combineObjFilesInSingleDir(QString qdir);
    static bool readWholeDir(const QString &qDirPath);
    static bool readSpecialDir(const QString &qDirPath);

    static bool test();
private:

};

#endif // OBJWORKER_H
