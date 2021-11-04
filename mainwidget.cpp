#include <mainwidget.h>

#include <QDebug>
#include <QHBoxLayout>

#include <iostream>
#include <fstream>

#include "glwidget.h"

#include "voxel_buffer.h"
#include "terraingenerator.h"
#include "meshGenerator/naiveSurfaceNets/surface_nets.h"

MainWidget::MainWidget(QWidget *parent) : QWidget(parent) {
    m_pGlWidget = new GlWidget(parent);
    m_pGlWidget->setMinimumSize(QSize(640, 480));

    QVBoxLayout *buttonLayout = new QVBoxLayout;
    buttonLayout->addWidget(new QLabel("x"));

    QHBoxLayout *mainlayout = new QHBoxLayout;
    mainlayout->addLayout(buttonLayout);
    mainlayout->addWidget(m_pGlWidget);

    setLayout(mainlayout);
}

MainWidget::~MainWidget() {}

void output2File(const TriMesh &mesh, std::string sOutputFile = "mesh.obj") {
    std::ofstream ofs(sOutputFile.c_str());
    char c[256];
    for (unsigned int xx = 0; xx < mesh.vertices.size(); ++xx) {
        sprintf(c, "v %f %f %f\n\0",
            mesh.vertices[xx].x,
            mesh.vertices[xx].y,
            mesh.vertices[xx].z);
        ofs << c;
        if (mesh.normals.size() > xx) {
            sprintf(c, "vn %f %f %f\n\0", mesh.normals[xx].x, mesh.normals[xx].y, mesh.normals[xx].z);
            ofs << c;
        }
    }
    for (unsigned int yy = 2; yy < mesh.indices.size(); yy += 3) {
        sprintf(c, "f %d %d %d\n\0", mesh.indices[yy - 2] + 1, mesh.indices[yy - 1] + 1, mesh.indices[yy] + 1);
        ofs << c;
    }
    ofs.close();
}

void test()
{
    std::shared_ptr<VoxelBuffer> buffer = std::make_shared<VoxelBuffer>();
    Vector3i vSize = Vector3i(100);
    buffer->create(100, 100, 100);

    TerrainGenerator generator;
    generator.generateTerrain(*buffer.get());

    auto const sdfFunction = [&](float x, float y, float z) -> float {
        return buffer->get_voxel_f(x, y, z, VoxelBuffer::CHANNEL_SDF);
    };

    SurfaceNets surfaceNets;
    TriMesh mesh = surfaceNets.surfaceNets(sdfFunction, vSize);

    output2File(mesh);
}
