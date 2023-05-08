#ifndef WORKTHREAD_H
#define WORKTHREAD_H

#include <QObject>
#include <QRunnable>

#include "meshGenerator/naiveSurfaceNets/voxel_mesher_surfacenet.h"

class WorkThread : public QObject, public QRunnable
{
    Q_OBJECT
public:
    struct InputBlock {
        std::shared_ptr<VoxelBuffer> voxels;
        Vector3i position;
        uint8_t lod = 0;
        bool can_be_discarded = true;
        float sort_heuristic = 0;
    };

    struct OutputBlock {
        VoxelMesher::Output blocky_surfaces;
        VoxelMesher::Output smooth_surfaces;
        Vector3i position;
        uint8_t lod = 0;
        bool drop_hint = false;
    };

    struct Input {
        std::vector<InputBlock> blocks;

        bool is_empty() const {
            return blocks.empty();
        }
    };

    struct Output {
        std::vector<OutputBlock> blocks;
    };

    void run() override {
        const InputBlock &block = input;
        OutputBlock &outputData = output;

        VoxelMesher::Input in = {*block.voxels, 0, input.position};
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
