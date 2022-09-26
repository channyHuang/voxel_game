#include <mainwidget.h>

#include <QDebug>
#include <QHBoxLayout>
#include <QPushButton>
#include <iostream>
#include <fstream>

#include "renderControl/glwidget.h"
#include "terrains/voxelterrain.h"
#include "commonSdf/sdf3.h"

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
        m_stBiomeParams.use_biomes = false;
        //TerrainGenerator_Roblox::getInstance()->generateTerrainByBiomes(tool, m_stBiomeParams);
    });

    QPushButton *testBtn = new QPushButton(this);
    testBtn->setText("only for test sub functions");
    connect(testBtn, &QPushButton::clicked, [&](){
        Sdf3 sdf3;
        sdf3.Sphere(10, Vector3(0, 0, 0));

        VoxelToolTerrain *tool = static_cast<VoxelToolTerrain*>(VoxelTerrain::getInstance()->get_voxel_tool());
        for (int i = -15; i <= 15; ++i) {
            for (int  j = -15; j <= 15; ++j) {
                for (int k = -15; k <= 15; ++k) {
                    tool->set_voxel_f(Vector3i(i, j, k), sdf3.getFun()(Vector3(i, j, k)));
                }
            }
        }

        VoxelMesherSurfaceNets surfaceNet;
        VoxelMesher::Output output;

        qDebug() << "testBtn end";
        //VoxelMesher::Input input({, 0, Vector3i(0, 0, 0)});
        //inpput.voxels;
        //surfaceNet.build(output, input);
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

