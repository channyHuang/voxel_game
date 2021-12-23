#include "objworker.h"

ObjWorker::ObjWorker()
{

}

ObjWorker::~ObjWorker()
{

}

bool ObjWorker::combineObjFilesInSingleDir(QString qdir) {
    QDir dir(qdir);
    if (!dir.exists()) return false;
    std::ofstream ofs(qdir.toStdString() + std::string(".obj"));
    qint64 offset = 0;
    QStringList qfilelist = dir.entryList(QDir::Files);
    foreach (QString qfile, qfilelist) {
        if (!qfile.endsWith(".obj")) continue;
        QFile file(dir.absolutePath() + "/" + qfile);
        if (!file.exists()) continue;
        if (!file.open(QFile::ReadOnly)) continue;

        qint64 count = 0;
        char line[256] = {0};
        qint64 linelength = file.readLine(line, 256);
        while (linelength != -1) {
            if (line[0] == 'v') {
                if (line[1] != 'n') {
                    count++;
                }
                ofs << line;
            } else { //f
                int a, b, c;
                QStringList faceinfo = QString(line).split(' ');
                a = faceinfo[1].toUInt();
                b = faceinfo[2].toUInt();
                c = faceinfo[3].toUInt();
                ofs << "f " << (a + offset) << " " << (b + offset) << " " << (c + offset) << std::endl;
            }

            linelength = file.readLine(line, 256);
        }
        file.close();
        offset += count;
    }
    ofs.close();
    return true;
}

bool ObjWorker::readWholeDir(const QString &qDirPath) {
    //read whole dir
    QDir wholedir(qDirPath);
    wholedir.setFilter(QDir::Dirs | QDir::NoDotAndDotDot);
    QFileInfoList dirlist = wholedir.entryInfoList();

    foreach(QFileInfo fileinfo, dirlist) {
        qDebug() << fileinfo.absolutePath() + "/" + fileinfo.fileName();
        if (!combineObjFilesInSingleDir(fileinfo.absolutePath() + "/" + fileinfo.fileName())) {
            return false;
        }
    }
    return true;
}

bool ObjWorker::readSpecialDir(const QString &qDirPath) {
    //read special dir
    std::vector<std::string> dirs = {"objFiles"};
    for (unsigned int i = 0; i < dirs.size(); i++) {
        QString qdirpath = qDirPath + "/" + QString(dirs[i].c_str());
        if (!combineObjFilesInSingleDir(qdirpath)) {
            return false;
        }
    }
    return true;
}

bool ObjWorker::test()
{
    QString sProPath = PRO_PATH;
    //test1
    readSpecialDir(sProPath);

    //test2
    //readWholeDir(sProPath + "/objFiles");

    qDebug() << "end";
    return true;
}
