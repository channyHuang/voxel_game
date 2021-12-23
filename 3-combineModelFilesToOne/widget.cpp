#include "widget.h"

#define DEFINE2STR(R) #R
#define STR(R) R

Widget::Widget(QWidget *parent)
    : QWidget(parent)
{
    comboBox = new QComboBox(this);
    comboBox->addItem(".obj");
    comboBox->addItem(".ply");

    pathLabel = new QLineEdit(this);

    pathBtn = new QPushButton(this);
    pathBtn->setText("select path");

    QHBoxLayout *pathLayout = new QHBoxLayout;
    pathLayout->addWidget(pathLabel);
    pathLayout->addWidget(pathBtn);

    QVBoxLayout *mainLayout = new QVBoxLayout;
    singleDirBtn = new QPushButton(this);
    singleDirBtn->setText("single path");

    wholeDirBtn = new QPushButton(this);
    wholeDirBtn->setText("whole path");

    mainLayout->addWidget(comboBox);
    mainLayout->addLayout(pathLayout);
    mainLayout->addWidget(singleDirBtn);
    mainLayout->addWidget(wholeDirBtn);

    setLayout(mainLayout);

    initConnection();
}

Widget::~Widget()
{
}

void Widget::initConnection() {
    connect(pathBtn, &QPushButton::clicked, this, &Widget::selectPath);
    connect(singleDirBtn, &QPushButton::clicked, this, &Widget::action);
    connect(wholeDirBtn, &QPushButton::clicked, [&](){
        QString dir = pathLabel->text();
        ObjWorker::readWholeDir(dir);
    });
}

void Widget::selectPath() {
    QString dir = QFileDialog::getExistingDirectory(this, tr("Open Directory"),
                                                    "C:",
                                                    QFileDialog::ShowDirsOnly
                                                    | QFileDialog::DontResolveSymlinks);
    pathLabel->setText(dir);
    update();
}

void Widget::action() {
    QString dir = pathLabel->text();
    switch(comboBox->currentIndex()) {
    case 0:
        ObjWorker::combineObjFilesInSingleDir(dir);
        break;
    case 1:
        PlyWorker::combinePlyFilesInSingleDir(dir);
        break;
    default:
        break;
    }
}
