#include "meshgeneratormanager.h"

MeshGeneratorManager::MeshGeneratorManager()
{
    pool.setMaxThreadCount(3);

    //sigMeshGenSuc.connect(this, &MeshGeneratorManager::sltFinish);
}

MeshGeneratorManager::~MeshGeneratorManager() {

}

void MeshGeneratorManager::sltMeshGenSuc(OutputBlock output) {
    std::unique_lock<std::mutex> lock(mutex);
    blocks.push_back(output);
    lock.unlock();
}

void MeshGeneratorManager::push(const Input& input) {
    for (int i = 0; i < input.blocks.size(); ++i) {
        OutputBlock output;

        const InputBlock &block = input.blocks[i];
        OutputBlock &outputData = output;

        MeshInput in = {*block.voxels, 0, block.position};
        output.position = block.position;
        VoxelMesherSurfaceNets* mesher = new VoxelMesherSurfaceNets;
        mesher->build(outputData.smooth_surfaces, in);

        sigMeshGenSuc(output);
    }
}

void MeshGeneratorManager::pop(Output &output) {
    output.blocks.swap(blocks);
    blocks.clear();
}
