#include <mainwidget.h>

#include <QDebug>
#include <QHBoxLayout>
#include <QPushButton>
#include <iostream>
#include <fstream>

#include "terrains/voxelBrush.h"
#include "terrains/terrainManager.h"
#include "voxelGenerator/terraingenerator_roblox.h"

MainWidget::MainWidget(QWidget *parent) : QWidget(parent) {
    init();

    initWidgets();
}

MainWidget::~MainWidget() {
    emit notice(TerrainEnum::Notification_Exit);
    if (terrain_thread != nullptr) {
        terrain_thread->exit();
    }
}

void MainWidget::initWidgets() {
    //m_pGlWidget = new GlWidget(this);
    //m_pGlWidget->setMinimumSize(QSize(640, 480));

    QPushButton *createBtn = new QPushButton(this);
    createBtn->setText("create");
    connect(createBtn, &QPushButton::clicked, [&](){
        VoxelBrush *brush = new VoxelBrush(TerrainManager::getInstance(), VoxelMap::getInstance());
        // generate terrain
        TerrainGenerator_Roblox::getInstance()->setRange(Vector3(0), m_size.to_vec3());
        TerrainGenerator_Roblox::getInstance()->generateTerrainByBiomes(brush, TerrainGenerator_Roblox::getInstance()->m_stBiomeParams);
    });

    //connect(TerrainManager::getInstance(), &TerrainManager::generateMeshSuc, m_pGlWidget, &GlWidget::updateMesh);

    // size
    QGroupBox *sizeBox = new QGroupBox("terrain size");
    QVBoxLayout *sizeLayout = new QVBoxLayout;
    for (int i = 0; i < 3; ++i) {
        m_pSizeBox[i] = new QSpinBox();
        m_pSizeBox[i]->setRange(0, 1024);
        m_pSizeBox[i]->setValue(10);

        connect(m_pSizeBox[i], &QSpinBox::valueChanged, [&, i](int value){
            m_size[i] = value;
        });
        sizeLayout->addWidget(m_pSizeBox[i]);
    }
    sizeBox->setLayout(sizeLayout);

    // biomes
    QGroupBox *biomeBox = new QGroupBox("terrain biomes");
    QVBoxLayout *biomeLayout = new QVBoxLayout;
    std::vector<std::string> biomeNames = {        "Water",
                                                   "Marsh",
                                                   "Plains",
                                                   "Hills",
                                                   "Dunes",
                                                   "Canyons",
                                                   "Mountains",
                                                   "Lavaflow",
                                                   "Arctic"};
    for (int i = 0; i < biomeNames.size(); ++i) {
        m_pBiomeBox[i] = new QCheckBox(biomeNames[i].data());

        connect(m_pBiomeBox[i], &QCheckBox::stateChanged, [&, i](int state){
            //if (state == Qt::Checked)
                //m_stBiomeParams.biomes_be_checked |= (1 << i);
            //else
                //m_stBiomeParams.biomes_be_checked ^= (1 << i);
        });
        biomeLayout->addWidget(m_pBiomeBox[i]);
    }
    biomeBox->setLayout(biomeLayout);

    QVBoxLayout *buttonLayout = new QVBoxLayout;
    buttonLayout->addWidget(sizeBox);
    buttonLayout->addWidget(biomeBox);
    buttonLayout->addWidget(createBtn);

    QHBoxLayout *mainlayout = new QHBoxLayout;
    mainlayout->addLayout(buttonLayout);
    //mainlayout->addWidget(m_pGlWidget);

    setLayout(mainlayout);
}

bool MainWidget::init() {
    TerrainManager *terrainManager = TerrainManager::getInstance();
    //connect(this, &MainWidget::notice, TerrainManager::getInstance(), &TerrainManager::_notification);
    //terrain_thread = new QThread(this);
    //terrain->moveToThread(terrain_thread);
    //terrainManager->moveToThread(terrain_thread);
    //terrain_thread->start(QThread::LowPriority);

    emit notice(TerrainEnum::Notification_Enter);
    timer.stop();
    timer.setInterval(10000);
    connect(&timer, &QTimer::timeout, [](){
        TerrainManager::getInstance()->_notification(TerrainEnum::Notification_Process);
    });

    timer.start(10000);
    return true;
}
