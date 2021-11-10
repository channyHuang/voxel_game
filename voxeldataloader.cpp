#include "voxeldataloader.h"

VoxelDataLoader::VoxelDataLoader() {
    // TODO I'm not sure it's worth to configure more than one thread for voxel streams

    FixedArray<Mgr::BlockProcessingFunc, Mgr::MAX_JOBS> processors;

    processors[0] = [this](ArraySlice<InputBlock> inputs, ArraySlice<OutputBlock> outputs, Mgr::ProcessorStats &stats) {
        this->process_blocks_thread_func(inputs, outputs, stats);
    };

    int batch_count = 128;
    int sync_interval_ms = 500;

    _mgr = new Mgr(1, sync_interval_ms, processors, true, batch_count);
}

VoxelDataLoader::~VoxelDataLoader() {
    if (_mgr) {
        delete(_mgr);
    }
}

// Can run in multiple threads
void VoxelDataLoader::process_blocks_thread_func(const ArraySlice<InputBlock> inputs, ArraySlice<OutputBlock> outputs, Mgr::ProcessorStats &stats) {
    std::vector<VoxelBlockRequest> emerge_requests;
    std::vector<VoxelBlockRequest> immerge_requests;

    for (size_t i = 0; i < inputs.size(); ++i) {

        const InputBlock &ib = inputs[i];

        int bs = 1 << _block_size_pow2;
        Vector3i block_origin_in_voxels = ib.position * (bs << ib.lod);

        if (ib.data.voxels_to_save == nullptr) {

            VoxelBlockRequest r;
            r.voxel_buffer->create(bs, bs, bs);
            r.origin_in_voxels = block_origin_in_voxels;
            r.lod = ib.lod;
            emerge_requests.push_back(r);

        } else {

            VoxelBlockRequest r;
            r.voxel_buffer = ib.data.voxels_to_save;
            r.origin_in_voxels = block_origin_in_voxels;
            r.lod = ib.lod;
            immerge_requests.push_back(r);
        }
    }

    // Assumes the stream won't change output order
    int iload = 0;
    for (size_t i = 0; i < outputs.size(); ++i) {

        const InputBlock &ib = inputs[i];
        OutputBlockData &output = outputs[i].data;

        if (ib.data.voxels_to_save == nullptr) {
            output.type = TYPE_LOAD;
            //output.voxels_loaded = emerge_requests.write[iload].voxel_buffer;
            ++iload;

        } else {
            output.type = TYPE_SAVE;
        }
    }

    // If unordered responses were allowed
    //
    //	size_t j = 0;
    //	for (size_t i = 0; i < emerge_requests.size(); ++i) {
    //		VoxelStream::BlockRequest &r = emerge_requests.write[i];
    //		OutputBlock &ob = outputs[j];
    //		ob.position = r.origin_in_voxels >> (_block_size_pow2 + r.lod);
    //		ob.lod = r.lod;
    //		ob.data.type = TYPE_LOAD;
    //		ob.data.voxels_loaded = r.voxel_buffer;
    //		++j;
    //	}
    //	for (size_t i = 0; i < immerge_requests.size(); ++i) {
    //		VoxelStream::BlockRequest &r = immerge_requests.write[i];
    //		OutputBlock &ob = outputs[j];
    //		ob.position = r.origin_in_voxels >> (_block_size_pow2 + r.lod);
    //		ob.lod = r.lod;
    //		ob.data.type = TYPE_SAVE;
    //		++j;
    //	}
}

