#ifndef VOXELDATALOADER_H
#define VOXELDATALOADER_H

#include "block_thread_manager.h"
#include "voxels/voxelBuffer.h"
#include "voxels/counting_time.h"

struct VoxelBlockRequest {
    VoxelBuffer* voxel_buffer;
    Vector3i origin_in_voxels;
    int lod;
};

class VoxelDataLoader {
public:
    struct InputBlockData {
        VoxelBuffer* voxels_to_save;
    };

    enum RequestType {
        TYPE_NOT_INITIALIZED = 0, // For error detection
        TYPE_SAVE,
        TYPE_LOAD
    };

    struct OutputBlockData {
        RequestType type;
        VoxelBuffer* voxels_loaded;
    };

    typedef VoxelBlockThreadManager<InputBlockData, OutputBlockData> Mgr;
    typedef Mgr::InputBlock InputBlock;
    typedef Mgr::OutputBlock OutputBlock;
    typedef Mgr::Input Input;
    typedef Mgr::Output Output;
    typedef Mgr::Stats Stats;

    VoxelDataLoader(unsigned int thread_count, unsigned int block_size_pow2);
    ~VoxelDataLoader();

    void push(const Input &input) { _mgr->push(input); }
    void pop(Output &output) { _mgr->pop(output); }

private:
    void process_blocks_thread_func(const ArraySlice<InputBlock> inputs, ArraySlice<OutputBlock> outputs, Mgr::ProcessorStats &stats);

    Mgr *_mgr = nullptr;
    int _block_size_pow2 = 0;
};


#endif // VOXELDATALOADER_H
