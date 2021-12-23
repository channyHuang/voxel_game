#ifndef HUMIDITY_H
#define HUMIDITY_H

#include "terrains/voxeltoolterrain.h"
#include "common_enum.h"

class Humidity
{
public:
    Humidity();
    ~Humidity();

    static Humidity* getInstance() {
        if (instance == nullptr) {
            instance = new Humidity();
        }
        return instance;
    }

    float getHumidity(VoxelToolTerrain* pVoxelTool, const Vector3i& vposi) {
        MaterialType material = (MaterialType)pVoxelTool->get_voxel(vposi, VoxelBuffer::CHANNEL_TYPE);
        switch(material) {
        case MaterialType::SNOW:
        case MaterialType::ICE:
            return 0.f;
        case MaterialType::WATER:
            return 0.2f;
        case MaterialType::SAND:
            return 1.f;
        default:
            break;
        }
        return 0.5f;
    }

private:
    static Humidity* instance;
};

#endif // HUMIDITY_H
