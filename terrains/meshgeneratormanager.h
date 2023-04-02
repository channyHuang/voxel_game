#ifndef MESHGENERATORMANAGER_H
#define MESHGENERATORMANAGER_H

#include <QObject>
#include <QThreadPool>
#include <QTimer>

#include <vector>
#include <unordered_map>
#include <algorithm>
#include <utility>

#include "workthread.h"

class MeshGeneratorManager : public QObject
{
    Q_OBJECT
public:
    explicit MeshGeneratorManager(QObject *parent = nullptr);
    ~MeshGeneratorManager();

    void push(const WorkThread::Input& input);
    void pop(WorkThread::Output &output);

    int get_minimum_padding() const { return _minimum_padding; }
    int get_maximum_padding() const { return _maximum_padding; }

signals:

public slots:
    void process();
    void sltFinish(WorkThread::OutputBlock output);

private:
    std::unordered_map<int, WorkThread> mapThreads;
    std::vector<WorkThread::OutputBlock> blocks;
    int index = 0;
    const int maxIndex = 10000;
    int _minimum_padding = 2;
    int _maximum_padding = 2;
    std::vector<WorkThread> threads;
    std::mutex mutex;
};

#endif // MESHGENERATORMANAGER_H
