#include "imageviewer.h"

ImageViewer::ImageViewer(QWidget *parent) : QWidget(parent)
{
    movie = new QMovie(this);
    imageLabel = new QLabel(this);

    selectBtn = new QPushButton(this);
    selectBtn->setText("select image");
    connect(selectBtn, &QPushButton::clicked, this, &ImageViewer::selectImage);

    fileLine = new QLineEdit(this);

    QHBoxLayout *fileLayout = new QHBoxLayout;
    fileLayout->addWidget(fileLine);
    fileLayout->addWidget(selectBtn);


    QVBoxLayout *mainLayout = new QVBoxLayout;
    mainLayout->addLayout(fileLayout);
    mainLayout->addWidget(imageLabel);

    setLayout(mainLayout);
}

void ImageViewer::selectImage() {
    QFileDialog dialog(this);
    QString qsFileName = dialog.getOpenFileName(this, tr(""), "", tr("Image files(*.png *.jpg *bmp *gif *ppm)"));

    if (qsFileName.endsWith(".gif")) {
        if (movie != nullptr) {
            movie->stop();
        }
        movie->setFileName(qsFileName);
        imageLabel->setMovie(movie);
        movie->start();
    } else {
        movie->stop();
        imageLabel->setPixmap(QPixmap::fromImage(QImage(qsFileName)));
    }
    update();
}

void ImageViewer::showEvent(QShowEvent *event) {
    if (movie && movie->isValid()) {
        movie->start();
    }
    QWidget::showEvent(event);
}

void ImageViewer::hideEvent(QHideEvent *event) {
    if (movie) {
        movie->stop();
    }
    QWidget::hideEvent(event);
}
