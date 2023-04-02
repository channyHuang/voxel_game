#ifndef WORKTHREAD_H
#define WORKTHREAD_H

#include <QObject>
#include <QRunnable>

#include "meshGenerator/naiveSurfaceNets/voxel_mesher_surfacenet.h"

class WorkThread : public QObject, public QRunnable
{
    Q_OBJECT
public:
    struct InputBlockData {
        VoxelBuffer* voxels;
    };

    struct OutputBlockData {
        VoxelMesher::Output blocky_surfaces;
        VoxelMesher::Output smooth_surfaces;
    };

    struct InputBlock {
        InputBlockData data;
        Vector3i position;
        uint8_t lod = 0;
        bool can_be_discarded = true;
        float sort_heuristic = 0;
    };

    struct OutputBlock {
        OutputBlockData data;
        Vector3i position;
        uint8_t lod = 0;
        bool drop_hint = false;
    };

    struct Input {
        std::vector<InputBlock> blocks;
        Vector3i priority_position;
        Vector3 priority_direction;
        int exclusive_region_extent = 0;
        int exclusive_region_max_lod = Math::MAX_LOD;
        bool use_exclusive_region = false;
        int max_lod_index = 0;

        bool is_empty() const {
            return blocks.empty();
        }
    };

    struct Output {
        std::vector<OutputBlock> blocks;
        // Stats stats;
    };



    void run() override {
        const InputBlockData &block = input.data;
        OutputBlockData &outputData = output.data;

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
