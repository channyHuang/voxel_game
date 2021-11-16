#include "voxel_mesh_updater.h"

VoxelMeshUpdater::VoxelMeshUpdater(unsigned int thread_count, MeshingParams params) {

    VoxelMesherSurfaceNets *mesher = new VoxelMesherSurfaceNets;

	_minimum_padding = 0;
	_maximum_padding = 0;


        _minimum_padding = max(_minimum_padding,
mesher->get_minimum_padding());
        _maximum_padding = max(_maximum_padding, mesher->get_maximum_padding());

	FixedArray<Mgr::BlockProcessingFunc, VoxelConstants::MAX_LOD> processors;

	for (unsigned int i = 0; i < thread_count; ++i) {

		if (i > 0) {

        }

        processors[i] = [this, mesher](const ArraySlice<InputBlock> inputs, ArraySlice<OutputBlock> outputs, Mgr::ProcessorStats &_) {
            this->process_blocks_thread_func(inputs, outputs, mesher);
		};
	}

    _mgr = new Mgr(thread_count, 50, processors);
}

VoxelMeshUpdater::~VoxelMeshUpdater() {
	if (_mgr) {
        delete(_mgr);
	}
}

void VoxelMeshUpdater::process_blocks_thread_func(
		const ArraySlice<InputBlock> inputs,
		ArraySlice<OutputBlock> outputs,
        VoxelMesher* mesher) {

	for (unsigned int i = 0; i < inputs.size(); ++i) {
		const InputBlock &ib = inputs[i];
		const InputBlockData &block = ib.data;
		OutputBlockData &output = outputs[i].data;

        VoxelMesher::Input input = {*block.voxels, ib.lod };
        input.position = ib.position;

        mesher->build(output.smooth_surfaces, input);
	}
}
