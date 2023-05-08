#ifndef WORKTHREAD_H
#define WORKTHREAD_H

#include <QObject>
#include <QRunnable>

#include "meshGenerator/naiveSurfaceNets/voxel_mesher_surfacenet.h"

#include "terrainCommonStruct.h"

class WorkThread : public QObject, public QRunnable
{
    Q_OBJECT
public:
    void run() override {
        const InputBlock &block = input;
        OutputBlock &outputData = output;

        MeshInput in = {*block.voxels, 0, input.position};
        output.position = input.position;
        mesher->build(outputData.smooth_surfaces, in);

        emit sigFinish(output);
    }
signals:
    void sigFinish(OutputBlock output);

public:
    int index;
    InputBlock input;
    OutputBlock output;
    VoxelMesherSurfaceNets *mesher = new VoxelMesherSurfaceNets;
};

#endif // WORKTHREAD_H
