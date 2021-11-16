#include <mainwidget.h>

#include <QDebug>
#include <QHBoxLayout>
#include <QPushButton>
#include <iostream>
#include <fstream>

#include "renderControl/glwidget.h"
#include "terrains/voxelterrain.h"

MainWidget::MainWidget(QWidget *parent) : QWidget(parent) {
    init();

    initWidgets();
}

MainWidget::~MainWidget() {
    emit notice(Notification_Exit);
    if (terrain_thread != nullptr) {
        terrain_thread->exit();
    }
}

void MainWidget::initWidgets() {
    m_pGlWidget = new GlWidget(this);
    m_pGlWidget->setMinimumSize(QSize(640, 480));

    QPushButton *createBtn = new QPushButton(this);
    createBtn->setText("create");
    connect(createBtn, &QPushButton::clicked, [&](){
        VoxelToolTerrain *tool = static_cast<VoxelToolTerrain*>(VoxelTerrain::getInstance()->get_voxel_tool());
        // generate terrain
        TerrainGenerator_Roblox::getInstance()->setRange(Vector3(0), Vector3(30, 30, 30));
        m_stBiomeParams.use_biomes = true;
        TerrainGenerator_Roblox::getInstance()->generateTerrainByBiomes(tool, m_stBiomeParams);
    });
    QVBoxLayout *buttonLayout = new QVBoxLayout;
    buttonLayout->addWidget(createBtn);

    QHBoxLayout *mainlayout = new QHBoxLayout;
    mainlayout->addLayout(buttonLayout);
    mainlayout->addWidget(m_pGlWidget);

    setLayout(mainlayout);
}

bool MainWidget::init() {
    auto terrain = VoxelTerrain::getInstance();
    connect(this, &MainWidget::notice, terrain, &VoxelTerrain::_notification);
    terrain_thread = new QThread(this);
    terrain->moveToThread(terrain_thread);
    terrain_thread->start(QThread::LowPriority);

    emit notice(Notification_Enter);
    timer.setInterval(5000);
    connect(&timer, &QTimer::timeout, terrain, [&](){
        terrain->_notification(Notification_Process);
    });
    timer.start();
    return true;
}

/*
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
*/

