#include "plyworker.h"

PlyWorker::PlyWorker()
{

}


PlyWorker::~PlyWorker()
{

}

bool PlyWorker::combinePlyFilesInSingleDir(QString qdir) {
    QDir dir(qdir);
    if (!dir.exists()) return false;
    std::ofstream ofs(qdir.toStdString() + std::string(".ply"));

    // calculate number of vertex and face
    qint64 totalVertex = 0, totalFace = 0, offset = 0;
    QStringList qfilelist = dir.entryList(QDir::Files);
    foreach (QString qfile, qfilelist) {
        if (!qfile.endsWith(".ply")) continue;
        QFile file(dir.absolutePath() + "/" + qfile);
        if (!file.exists()) continue;
        if (!file.open(QFile::ReadOnly | QIODevice::Text)) continue;

        QTextStream in(&file);
        QString line = in.readLine();
        while (!line.isNull()) {
            if (line.startsWith("end_header")) {
                break;
            }
            if (line.startsWith("element vertex")) {
                totalVertex += line.remove("element vertex ").toUInt();
            } else if (line.startsWith("element face")) {
                totalFace += line.remove("element face ").toUInt();
            }
            line = in.readLine();
        }
        file.close();
    }
    qDebug() << "total " << totalVertex << " " << totalFace;

    // write header
    ofs << "ply\nformat ascii 1.0\n";
    ofs << "element vertex " << totalVertex << std::endl;
    ofs << "property float x\nproperty float y\nproperty float z\nproperty uchar red\nproperty uchar green\nproperty uchar blue\nproperty float nx\nproperty float ny\nproperty float nz\n";
    ofs << "element face " << totalFace << std::endl;
    ofs << "property list uchar int vertex_index\nend_header\n";

    // write vertex
    foreach (QString qfile, qfilelist) {
        if (!qfile.endsWith(".ply")) continue;
        QFile file(dir.absolutePath() + "/" + qfile);
        if (!file.exists()) continue;
        if (!file.open(QFile::ReadOnly | QIODevice::Text)) continue;

        QTextStream in(&file);
        QString line = in.readLine();
        qint64 vertex = 0;
        while (!line.isNull()) {
            if (line.startsWith("end_header")) {
                break;
            }
            if (line.startsWith("element vertex")) {
                vertex = line.remove("element vertex ").toUInt();
            }
            line = in.readLine();
        }

        if (vertex == 0) continue;
        line = in.readLine();
        while (!line.isNull()) {
            if (vertex <= 0) {
                break;
            }
            vertex--;
            ofs << line.toStdString() << "\n";
            line = in.readLine();
        }
        file.close();
    }

    // write face
    foreach (QString qfile, qfilelist) {
        if (!qfile.endsWith(".ply")) continue;
        QFile file(dir.absolutePath() + "/" + qfile);
        if (!file.exists()) continue;
        if (!file.open(QFile::ReadOnly | QIODevice::Text)) continue;

        QTextStream in(&file);
        QString line = in.readLine();
        qint64 vertex = 0, face = 0;
        while (!line.isNull()) {
            if (line.startsWith("end_header")) {
                break;
            }
            if (line.startsWith("element vertex")) {
                vertex = line.remove("element vertex ").toUInt();
            } else if (line.startsWith("element face")) {
                face = line.remove("element face ").toUInt();
            }
            line = in.readLine();
        }

        if (face == 0) continue;
        line = in.readLine();
        qint64 curvertex = vertex;
        while (!line.isNull()) {
            if (curvertex > 0) {
                curvertex--;
                line = in.readLine();
                continue;
            }
            if (face <= 0) break;
            face--;
            QStringList qslIndex = line.split(' ');
            ofs << qslIndex[0].toStdString();
            for (int i = 1; i < qslIndex.size(); ++i) {
                qint64 index = qslIndex[i].toUInt();
                ofs << " " << index + offset;
            }
            ofs << "\n";

            line = in.readLine();
        }

        offset += vertex;
        file.close();
    }

    ofs.close();
    return true;
}

bool PlyWorker::readWholeDir(const QString &qDirPath) {
    //read whole dir
    QDir wholedir(qDirPath);
    wholedir.setFilter(QDir::Dirs | QDir::NoDotAndDotDot);
    QFileInfoList dirlist = wholedir.entryInfoList();

    foreach(QFileInfo fileinfo, dirlist) {
        qDebug() << fileinfo.absolutePath() + "/" + fileinfo.fileName();
        if (!combinePlyFilesInSingleDir(fileinfo.absolutePath() + "/" + fileinfo.fileName())) {
            return false;
        }
    }
    return true;
}

bool PlyWorker::readSpecialDir(const QString &qDirPath) {
    //read special dir
    std::vector<std::string> dirs = {"objFiles"};
    for (unsigned int i = 0; i < dirs.size(); i++) {
        QString qdirpath = qDirPath + "/" + QString(dirs[i].c_str());
        if (!combinePlyFilesInSingleDir(qdirpath)) {
            return false;
        }
    }
    return true;
}

bool PlyWorker::test()
{
    QString sProPath = PRO_PATH;
    //test1
    readSpecialDir(sProPath);

    //test2
    //readWholeDir(sProPath + "/objFiles");

    qDebug() << "end";
    return true;
}
