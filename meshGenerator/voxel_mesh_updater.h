#ifndef VOXEL_MESH_UPDATER_H
#define VOXEL_MESH_UPDATER_H

#include "../voxels/voxel_buffer.h"
#include "../meshGenerator/voxel_mesher.h"

#include "terrains/block_thread_manager.h"
#include "meshGenerator/voxel_mesher.h"
#include "meshGenerator/naiveSurfaceNets/voxel_mesher_surfacenet.h"

class VoxelMeshUpdater {
public:
    struct InputBlockData {
        VoxelBuffer* voxels;
    };

    struct OutputBlockData {
        VoxelMesher::Output blocky_surfaces;
        VoxelMesher::Output smooth_surfaces;
    };

    struct MeshingParams {
        //Ref<VoxelLibrary> library;
        bool baked_ao = true;
        float baked_ao_darkness = 0.8;
        bool smooth_surface = false;
    };

    typedef VoxelBlockThreadManager<InputBlockData, OutputBlockData> Mgr;
    typedef Mgr::InputBlock InputBlock;
    typedef Mgr::OutputBlock OutputBlock;
    typedef Mgr::Input Input;
    typedef Mgr::Output Output;
    typedef Mgr::Stats Stats;

    VoxelMeshUpdater(unsigned int thread_count, MeshingParams params);
    ~VoxelMeshUpdater();

    void push(const Input &input) { _mgr->push(input); }
    void pop(Output &output) { _mgr->pop(output); }

    int get_minimum_padding() const { return _minimum_padding; }
    int get_maximum_padding() const { return _maximum_padding; }

private:
    void process_blocks_thread_func(const ArraySlice<InputBlock> inputs,
            ArraySlice<OutputBlock> outputs,
            VoxelMesher *mesher);

    Mgr *_mgr = nullptr;
    int _minimum_padding = 0;
    int _maximum_padding = 0;
};

#endif // VOXEL_MESH_UPDATER_H
