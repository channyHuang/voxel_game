#ifndef PLYWORKER_H
#define PLYWORKER_H

#include <vector>
#include <QDir>
#include <QString>
#include <QDebug>
#include <QFile>
#include <QFileInfo>
#include <fstream>

class PlyWorker
{
public:
    PlyWorker();
    ~PlyWorker();

    static bool combinePlyFilesInSingleDir(QString qdir);
    static bool readWholeDir(const QString &qDirPath);
    static bool readSpecialDir(const QString &qDirPath);

    static bool test();
};

#endif // PLYWORKER_H
