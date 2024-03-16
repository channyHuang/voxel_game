#ifndef TERRAINCOMMONSTRUCT_H
#define TERRAINCOMMONSTRUCT_H

#include <vector>
#include <array>

#include "voxels/voxelBuffer.h"

namespace TerrainStruct {

struct Arrays {
    std::vector<Vector3> positions;
    std::vector<Vector3> normals;
    std::vector<uint32_t> indices;
    std::vector<uint32_t> materials;
    std::array<std::vector<uint32_t>, /*Cube::SIDE_COUNT*/6> transition_surfaces;
    std::vector<uint32_t> indices_water;
    bool isWater = false;

    bool empty() {
        return (positions.size() <= 0 || indices.size() <= 0);
    }
};

struct MeshInput {
    const VoxelBuffer &voxels;
    int lod;
    Vector3i position;
};

struct MeshOutput {
    std::vector<Arrays> surfaces;
};

struct InputBlock {
    std::shared_ptr<VoxelBuffer> voxels;
    Vector3i position;
    uint8_t lod = 0;
    bool can_be_discarded = true;
    float sort_heuristic = 0;
};

struct OutputBlock {
    MeshOutput blocky_surfaces;
    MeshOutput smooth_surfaces;
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

}

#endif // TERRAINCOMMONSTRUCT_H
