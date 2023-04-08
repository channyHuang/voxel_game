#ifndef WIDGET_H
#define WIDGET_H

#include <QOpenGLWidget>
#include <QOpenGLBuffer>
#include <QOpenGLFunctions>
#include <QOpenGLShaderProgram>
#include <QOpenGLVertexArrayObject>
#include <QOpenGLTexture>
#include <QMouseEvent>
#include <QMatrix4x4>
#include <QTimer>

#include <unordered_map>

#include "meshGenerator/voxel_mesher.h"
#include "objfileloader.h"
#include "vector2.h"
#include "common_enum.h"

class GlWidget : public QOpenGLWidget
{
    Q_OBJECT
public:
    struct MeshData {
        Vector3 positions;
        Vector3 normals;
        Vector2 uvs;
        uint32_t materials = 1;

        MeshData(const Vector3& pos, const Vector3& n, const Vector2& uv, uint32_t mat = 0)
            : positions(pos), normals(n), uvs(uv), materials(mat) {}
    };

    GlWidget(QWidget *partent = nullptr);
    ~GlWidget();

public slots:
    void updateMesh(const Arrays& surface, const Vector3i& pos);
protected:
    void initializeGL() override;
    void paintGL() override;
    void resizeGL(int w, int h) override;

    void mousePressEvent(QMouseEvent *e) override;
    void mouseMoveEvent(QMouseEvent *e) override;
    void mouseReleaseEvent(QMouseEvent *e) override;
    void keyPressEvent(QKeyEvent *e) override;
    void keyReleaseEvent(QKeyEvent *e) override;
    void wheelEvent(QWheelEvent *e) override;
private:
    void initTexture();
    void drawCubeGeometry();
    void initShader();
    void setRotation(int angle, int axis);
    void bindTexture();

private:
    QString sProPath;
    QOpenGLShaderProgram *m_shader;
    QOpenGLVertexArrayObject *m_vao;
    QOpenGLFunctions *f = nullptr;
    QOpenGLTexture *texture;
    std::unordered_map<MaterialType, QOpenGLTexture *> textures;

    int m_projMatrixLoc;
    int m_mvMatrixLoc;
    int m_normalMatrixLoc;
    int m_lightPosLoc;

    QMatrix4x4 projection;
    QVector2D mousePressPosition;

    qreal zNear = 0.001f, zFar = 1000.f, fov = 60.f, aspect;
    bool m_core = false;

    ObjFileLoader m_dataLoader;
    QOpenGLBuffer *m_vbo;

    int m_Rot[3] = {0, 0, 0};
    QMatrix4x4 m_camera;
    QMatrix4x4 m_world;
    Vector3 m_vCameraDir = Vector3(0, 0, 1);
    Vector3 m_vCameraUp = Vector3(0, 1, 0);
    Vector3 m_vCameraPos = Vector3(0, 50, -50);

    QTimer timer;
    int m_nNumOfTris = 0;
    std::unordered_map<Vector3i, std::vector<MeshData>, Vector3iHash> mapMeshs;
    std::mutex mutex;
};

#endif
