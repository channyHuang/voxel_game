#ifndef MESHGENERATORMANAGER_H
#define MESHGENERATORMANAGER_H

#include <vector>
#include <unordered_map>
#include <algorithm>
#include <utility>

#include "commonFunc/threadPool.h"
#include "commonFunc/signalSlots.h"

#include "terrainCommonStruct.h"
using namespace TerrainStruct;

class MeshGeneratorManager
{
public:
    explicit MeshGeneratorManager();
    ~MeshGeneratorManager();

    void push(const Input& input);
    void pop(Output &output);

    int get_minimum_padding() const { return _minimum_padding; }
    int get_maximum_padding() const { return _maximum_padding; }

// signals:
    SignalSlots::Signal<void(OutputBlock)> sigMeshGenSuc;

// public slots:
    void sltMeshGenSuc(OutputBlock output);

private:
    std::vector<OutputBlock> blocks;
    int index = 0;
    const int maxIndex = 10000;
    int _minimum_padding = 2;
    int _maximum_padding = 2;
    std::mutex mutex;

    ThreadPool pool;
};

#endif // MESHGENERATORMANAGER_H
