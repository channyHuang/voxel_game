#include "meshgeneratormanager.h"

MeshGeneratorManager::MeshGeneratorManager(QObject *parent)
    : QObject{parent}
{
    QThreadPool::globalInstance()->setMaxThreadCount(3);
    qRegisterMetaType<WorkThread::OutputBlock>("OutputBlock");
    qRegisterMetaType<Arrays>("Arrays");
    qRegisterMetaType<Vector3i>("Vector3i");
}

MeshGeneratorManager::~MeshGeneratorManager() {
    QThreadPool::globalInstance()->waitForDone();
}

void MeshGeneratorManager::process() {

}

void MeshGeneratorManager::sltFinish(WorkThread::OutputBlock output) {
    std::unique_lock<std::mutex> lock(mutex);
    blocks.push_back(output);
    lock.unlock();
}

void MeshGeneratorManager::push(const WorkThread::Input& input) {
    std::vector<WorkThread> tmpThread(input.blocks.size());
    threads.swap(tmpThread);
    for (int i = 0; i < input.blocks.size(); ++i) {

        threads[i].index = index++;
        if (index >= maxIndex) {
            index = 0;
        }
        threads[i].setAutoDelete(false);
        threads[i].input = input.blocks[i];
        connect(&threads[i], &WorkThread::sigFinish, this, &MeshGeneratorManager::sltFinish);
        QThreadPool::globalInstance()->start(&threads[i]);
    }
}

void MeshGeneratorManager::pop(WorkThread::Output &output) {
    output.blocks.swap(blocks);
    blocks.clear();
}
