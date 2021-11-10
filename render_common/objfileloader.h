#ifndef OBJFILELOADER_H
#define OBJFILELOADER_H

#include <iostream>
#include <fstream>
#include <vector>
#include <QVector3D>

#include <QOpenGLFunctions>

struct FileObject {
    std::vector<QVector3D> positions;
    std::vector<QVector3D> faces;
    std::vector<QVector3D> normals;
    std::vector<QVector3D> face_normals;

    void clear() {
        positions.clear();
        faces.clear();
        face_normals.clear();
        normals.clear();
    }
};

class ObjFileLoader
{
public:
    ObjFileLoader();

    const GLfloat *constData() const { return m_data.constData(); }
    int count() const { return m_count; }
    int vertexCount() const { return m_count / 6; }

    bool loadObjFile(std::string sObjFile);
    bool loadObjFile(std::string sObjFile, FileObject &stFileObject);
    bool pushToData(FileObject &stFileObject);

private:
    QVector<GLfloat> m_data;
    int m_count = 0;

    FileObject m_stFileObject;
};

#endif // OBJFILELOADER_H
