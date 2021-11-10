#pragma once

#include "FastNoiseLite.h"
#include "voxel_buffer.h"
#include "common_math/box.h"

class TerrainGenerator
{
public:
    enum TerrainType {
        Terrain_Default = 0,
        Terrain_Flat,
        Terrain_Moutain
    };

    TerrainGenerator();
    ~TerrainGenerator();

    void generateTerrain(VoxelBuffer &buffer, TerrainType eTerrainType = Terrain_Default);

private:
    Box bRange;
    FastNoiseLite noise_lite;
};
