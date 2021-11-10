#ifndef FLIPVIEW_H
#define FLIPVIEW_H

#include <QWidget>
#include <QHBoxLayout>
#include <QVBoxLayout>
#include <QLabel>
#include <QPixmap>
#include <QImage>

class FlipView : public QWidget
{
    Q_OBJECT
public:
    explicit FlipView(QWidget *parent = nullptr);

signals:
private:
    QLabel *background_label = nullptr;
};

#endif // FLIPVIEW_H
