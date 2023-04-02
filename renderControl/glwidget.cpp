#include "glwidget.h"

#define CONFIG2STR(R) #R
#define CONFIG2QSTR(R) CONFIG2STR(R)

#include "uvUnwrap/xatlas.h"

GlWidget::GlWidget(QWidget* parent):
    QOpenGLWidget (parent),
    m_shader(nullptr), m_vao(nullptr), m_vbo(nullptr)
{
    m_core = QSurfaceFormat::defaultFormat().profile() == QSurfaceFormat::CoreProfile;
    sProPath = CONFIG2QSTR(PRO_PATH);
    setFocusPolicy(Qt::ClickFocus);
}

GlWidget::~GlWidget(){

}

static void qNormalizeAngle(int &angle)
{
    while (angle < 0)
        angle += 360 * 16;
    while (angle > 360 * 16)
        angle -= 360 * 16;
}

void GlWidget::setRotation(int angle, int axis)
{
    qNormalizeAngle(angle);
    if (angle != m_Rot[axis]) {
        m_Rot[axis] = angle;
        update();
    }
}

void GlWidget::mousePressEvent(QMouseEvent *e) {
    mousePressPosition = QVector2D(e->localPos());
    QOpenGLWidget::mousePressEvent(e);
}

void GlWidget::mouseMoveEvent(QMouseEvent *e) {
    int dx = e->x() - mousePressPosition.x();
    int dy = e->y() - mousePressPosition.y();

    if (e->buttons() & Qt::LeftButton) {
        setRotation(m_Rot[0] + 8 * dy, 0);
        setRotation(m_Rot[1] + 8 * dx, 1);
    } else if (e->buttons() & Qt::RightButton) {
        setRotation(m_Rot[0] + 8 * dy, 0);
        setRotation(m_Rot[2] + 8 * dx, 2);
    }
    mousePressPosition = QVector2D(e->localPos());
    update();

    QOpenGLWidget::mouseMoveEvent(e);
}

void GlWidget::mouseReleaseEvent(QMouseEvent *e) {
    QOpenGLWidget::mouseReleaseEvent(e);
}

void GlWidget::keyReleaseEvent(QKeyEvent *e) {
    switch(e->key()) {
    case Qt::Key_W:
        m_vCameraPos += m_vCameraDir;
        break;
    case Qt::Key_S:
        m_vCameraPos -= m_vCameraDir;
        break;
    case Qt::Key_A:
        m_vCameraPos += m_vCameraUp.cross(m_vCameraDir);
        break;
    case Qt::Key_D:
        m_vCameraPos -= m_vCameraUp.cross(m_vCameraDir);
        break;
    case Qt::Key_Q:
        m_vCameraPos += m_vCameraUp;
        break;
    case Qt::Key_E:
        m_vCameraPos -= m_vCameraUp;
        break;
    default:
        break;
    }

    m_camera.setToIdentity();
    m_camera.translate(m_vCameraPos.x, m_vCameraPos.y, m_vCameraPos.z);
    update();

    QOpenGLWidget::keyReleaseEvent(e);
}

void GlWidget::wheelEvent(QWheelEvent *e) {
    m_camera.setToIdentity();

    if (e->angleDelta().y() > 0) {
        m_vCameraPos += m_vCameraDir;
    } else {
        m_vCameraPos -= m_vCameraDir;
    }

    m_camera.translate(m_vCameraPos.x, m_vCameraPos.y, m_vCameraPos.z);
    update();
}

void GlWidget::initShader() {
    m_shader = new QOpenGLShaderProgram;
    m_shader->addShaderFromSourceFile(QOpenGLShader::Vertex, m_core ? ":/vshader_core.glsl" : ":/vshader.glsl");
    m_shader->addShaderFromSourceFile(QOpenGLShader::Fragment, m_core ? ":/fshader_core.glsl" : ":/fshader.glsl");

    m_shader->bindAttributeLocation("vertex", 0);
    m_shader->bindAttributeLocation("normal", 1);
    m_shader->bindAttributeLocation("a_texcoord", 2);

    if (m_shader->link()) {
        qDebug("Shaders link success.");
    } else {
        qDebug("Shaders link failed!");
    }

    m_shader->bind();
    m_projMatrixLoc = m_shader->uniformLocation("projMatrix");
    m_mvMatrixLoc = m_shader->uniformLocation("mvMatrix");
    m_normalMatrixLoc = m_shader->uniformLocation("normalMatrix");
    m_lightPosLoc = m_shader->uniformLocation("lightPos");
}

void GlWidget::initializeGL()
{
    //init gl environment
    f = this->context()->functions();
    f->initializeOpenGLFunctions();
    f->glClearColor(0.0f, 0.0f, 0.0f, 1.0f);

    initShader();
    initTexture();

    m_vao = new QOpenGLVertexArrayObject();
    m_vbo = new QOpenGLBuffer(QOpenGLBuffer::Type::VertexBuffer);
    m_vao->create();
    m_vbo->create();

    m_vao->bind();
    m_vbo->bind();
    //m_dataLoader.loadObjFile((sProPath + "/gazebo.obj").toStdString());
    //m_vbo->allocate(m_dataLoader.constData(), m_dataLoader.count() * sizeof(GLfloat));

    f->glEnableVertexAttribArray(0);
    f->glEnableVertexAttribArray(1);
    f->glEnableVertexAttribArray(2);
    f->glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 8 * sizeof(GLfloat), 0);
    f->glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 8 * sizeof(GLfloat), reinterpret_cast<void *>(3 * sizeof(GLfloat)));
    f->glVertexAttribPointer(2, 2, GL_FLOAT, GL_FALSE, 8 * sizeof(GLfloat), reinterpret_cast<void *>(6 * sizeof(GLfloat)));

    m_camera.setToIdentity();
    m_camera.translate(0, 0, -1);
    m_shader->setUniformValue(m_lightPosLoc, QVector3D(0, 0, 70));

    m_shader->release();
    m_vbo->release();
    m_vao->release();
}

void GlWidget::initTexture() {
    texture = new QOpenGLTexture(QImage(":/cork.jpg").mirrored());
    texture->setMinificationFilter(QOpenGLTexture::Nearest);
    texture->setMagnificationFilter(QOpenGLTexture::Linear);
    texture->setWrapMode(QOpenGLTexture::Repeat);


}

void GlWidget::paintGL()
{
    //qDebug() << __FUNCTION__;
    QOpenGLFunctions *f = this->context()->functions();
    f->glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    f->glEnable(GL_DEPTH_TEST);
    f->glEnable(GL_CULL_FACE);

    m_world.setToIdentity();
    m_world.rotate(180.0f - (m_Rot[0] / 16.0f), 1, 0, 0);
    m_world.rotate(m_Rot[1] / 16.0f, 0, 1, 0);
    m_world.rotate(m_Rot[2] / 16.0f, 0, 0, 1);

    m_vao->bind();
    m_vbo->bind();
    m_shader->bind();

    std::unique_lock<std::mutex> lock(mutex);
    QVector<GLfloat> data(m_nNumOfTris * 8);
    GLfloat *p = data.data();
    int count = m_nNumOfTris;
    for (auto itr = mapMeshs.begin(); itr != mapMeshs.end(); itr++) {
        for (int j = 0; j < itr->second.size(); j++) {
            *p++ = itr->second[j].positions.x;
            *p++ = itr->second[j].positions.y;
            *p++ = itr->second[j].positions.z;
            *p++ = itr->second[j].normals.x;
            *p++ = itr->second[j].normals.y;
            *p++ = itr->second[j].normals.z;
            *p++ = itr->second[j].uvs.x;
            *p++ = itr->second[j].uvs.y;
        }
    }
    lock.unlock();
    m_vbo->allocate(data.constData(), m_nNumOfTris * 8 * sizeof(GLfloat));

    texture->bind();
    m_shader->setUniformValue("texture", 0);

    m_shader->setUniformValue(m_projMatrixLoc, projection);
    m_shader->setUniformValue(m_mvMatrixLoc, m_camera * m_world);
    QMatrix3x3 normalMatrix = m_world.normalMatrix();
    m_shader->setUniformValue(m_normalMatrixLoc, normalMatrix);

    f->glEnableVertexAttribArray(0);
    f->glEnableVertexAttribArray(1);
    f->glEnableVertexAttribArray(2);
    f->glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 8 * sizeof(GLfloat), 0);
    f->glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 8 * sizeof(GLfloat), reinterpret_cast<void *>(3 * sizeof(GLfloat)));
    f->glVertexAttribPointer(2, 2, GL_FLOAT, GL_FALSE, 8 * sizeof(GLfloat), reinterpret_cast<void *>(6 * sizeof(GLfloat)));


    f->glDrawArrays(GL_TRIANGLES, 0, count);

//    f->glDrawArrays(GL_TRIANGLES, 0, m_dataLoader.vertexCount());
    m_shader->release();
    m_vbo->release();
    m_vao->release();
}

void GlWidget::resizeGL(int w, int h)
{
    aspect = qreal(w) / qreal(h ? h : 1);
    projection.setToIdentity();
    projection.perspective(fov, aspect, zNear, zFar);
}

void output(const Arrays& surface, const Vector3i& position) {
    std::string sFileName = "glwidget_" + std::to_string(position.x) + "_" + std::to_string(position.y) + "_" + std::to_string(position.z) + ".obj";
    std::ofstream ofs(sFileName);
    for (size_t i = 0; i < surface.positions.size(); i++) {
        ofs << "v " << surface.positions[i][0] << " " << surface.positions[i][1]  << " " << surface.positions[i][2] << std::endl;
        ofs << "vn " << surface.normals[i].x << " " << surface.normals[i].y << " " << surface.normals[i].z << std::endl;
    }
    for (size_t i = 0; i < surface.indices.size(); i += 3) {
        ofs << "f " << surface.indices[i] + 1 << " " << surface.indices[i + 2] + 1 << " " << surface.indices[i + 1] + 1 << std::endl;
    }
    ofs.close();
}

void GlWidget::updateMesh(Arrays surface, Vector3i pos) {
    //qDebug() << __FUNCTION__ << " " << pos.toString().c_str();
    xatlas::Atlas *atlas = xatlas::Create();
    xatlas::MeshDecl meshDecl;
    meshDecl.vertexCount = surface.positions.size();
    meshDecl.vertexPositionData = surface.positions.data();
    meshDecl.vertexPositionStride = sizeof(float) * 3;
    meshDecl.vertexNormalData = surface.normals.data();
    meshDecl.vertexNormalStride = sizeof(float) * 3;
    meshDecl.indexCount = surface.indices.size();
    meshDecl.indexData = surface.indices.data();
    meshDecl.indexFormat = xatlas::IndexFormat::UInt32;
    xatlas::AddMeshError error = xatlas::AddMesh(atlas, meshDecl, (uint32_t)1);
    if (error != xatlas::AddMeshError::Success) {
        qDebug() << __FUNCTION__ << " uvUnwrap failed";
        xatlas::Destroy(atlas);
        return;
    }
    xatlas::Generate(atlas);

    std::unique_lock<std::mutex> lock(mutex);
    auto itr = mapMeshs.find(pos);
    if (itr == mapMeshs.end()) {
        mapMeshs[pos] = std::vector<MeshData>();
    } else {
        m_nNumOfTris -= mapMeshs[pos].size();
        mapMeshs[pos].clear();
    }

    if (atlas->meshCount <= 0) return;
    const xatlas::Mesh &mesh = atlas->meshes[0];
    for (int i = 0; i < surface.indices.size(); i += 3) {
        xatlas::Vertex &vertex = mesh.vertexArray[surface.indices[i]];

        mapMeshs[pos].push_back(MeshData(surface.positions[surface.indices[i]],
                                surface.normals[surface.indices[i]],
                                Vector2(vertex.uv[0], vertex.uv[1])));

        vertex = mesh.vertexArray[surface.indices[i + 2]];
        mapMeshs[pos].push_back(MeshData(surface.positions[surface.indices[i + 2]],
                                surface.normals[surface.indices[i + 2]],
                                Vector2(vertex.uv[0], vertex.uv[1])));

        vertex = mesh.vertexArray[surface.indices[i + 1]];
        mapMeshs[pos].push_back(MeshData(surface.positions[surface.indices[i + 1]],
                                surface.normals[surface.indices[i + 1]],
                                Vector2(vertex.uv[0], vertex.uv[1])));
    }
    m_nNumOfTris += surface.indices.size();
    lock.unlock();
    //output(surface, pos);

    xatlas::Destroy(atlas);
}

