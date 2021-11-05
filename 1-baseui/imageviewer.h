#ifndef IMAGEVIEWER_H
#define IMAGEVIEWER_H

#include <QWidget>
#include <QVBoxLayout>
#include <QLabel>
#include <QPushButton>
#include <QFileDialog>
#include <QLineEdit>
#include <QPixmap>
#include <QMovie>

class ImageViewer : public QWidget
{
    Q_OBJECT
public:
    explicit ImageViewer(QWidget *parent = nullptr);

signals:

public slots:
    void selectImage();

protected:
    void showEvent(QShowEvent *event) override;
    void hideEvent(QHideEvent *event) override;

private:
    QLabel *imageLabel = nullptr;
    QPushButton *selectBtn = nullptr;
    QLineEdit *fileLine = nullptr;
    QMovie *movie = nullptr;
};

#endif // IMAGEVIEWER_H
