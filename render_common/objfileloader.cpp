#include "objfileloader.h"

ObjFileLoader::ObjFileLoader()
{
    m_stFileObject.clear();
}

bool ObjFileLoader::loadObjFile(std::string sObjFile) {
    bool res = loadObjFile(sObjFile, m_stFileObject);
    pushToData(m_stFileObject);
    return res;
}

bool ObjFileLoader::loadObjFile(std::string sObjFile, FileObject &stFileObject) {
    std::ifstream ifs(sObjFile);
    if (!ifs.is_open()) {
        return false;
    }
    int count = stFileObject.positions.size();
    QVector3D vec, veci[3];
    std::string value;
    char tmp[10];
    while (std::getline(ifs, value)) {
        if (value.length() <= 5) continue;
        switch (value[0]) {
        case 'f':
            sscanf(value.c_str(), "%s %f/%f/%f %f/%f/%f %f/%f/%f", &tmp,
                   &veci[0][0], &veci[1][0], &veci[2][0],
                    &veci[0][1], &veci[1][1], &veci[2][1],
                    &veci[0][2], &veci[1][2], &veci[2][2]);

            stFileObject.faces.push_back(veci[0] - QVector3D(1, 1, 1));
            stFileObject.face_normals.push_back(veci[2] - QVector3D(1, 1, 1));
            break;
        case 'v':
            sscanf(value.c_str(), "%s %f %f %f", tmp, &vec[0], &vec[1], &vec[2]);
            if (value[1] == 'n') {
                stFileObject.normals.push_back(vec);
            } else if (value[1] == ' ') {
                stFileObject.positions.push_back(vec + QVector3D(count, count, count));
            }
            break;
        default:
            break;
        }
    }
    ifs.close();
    return (stFileObject.positions.size() == stFileObject.normals.size());
}

bool ObjFileLoader::pushToData(FileObject &stFileObject) {
    m_data.resize(stFileObject.faces.size() * 3 * 6);
    m_count = 0;

    for (size_t i = 0; i < stFileObject.faces.size(); ++i) {
        GLfloat *p = m_data.data() + m_count;
        *p++ = stFileObject.positions[stFileObject.faces[i][0]].x();
        *p++ = stFileObject.positions[stFileObject.faces[i][1]].y();
        *p++ = stFileObject.positions[stFileObject.faces[i][2]].z();
        *p++ = stFileObject.normals[stFileObject.face_normals[i][0]].x();
        *p++ = stFileObject.normals[stFileObject.face_normals[i][1]].y();
        *p++ = stFileObject.normals[stFileObject.face_normals[i][2]].z();
        m_count += 6;
    }
    return true;
}
