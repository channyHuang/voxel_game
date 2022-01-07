#ifndef TERRAINGENERATOR_H
#define TERRAINGENERATOR_H

#include "FastNoiseLite.h"
#include "voxel_buffer.h"
#include "box.h"

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

#endif
