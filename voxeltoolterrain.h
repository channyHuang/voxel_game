#ifndef VOXELTOOLTERRAIN_H
#define VOXELTOOLTERRAIN_H


#include "voxeltool.h"
#include "voxelmap.h"
class VoxelTerrain;

class VoxelToolTerrain : public VoxelTool {
public:
    VoxelToolTerrain();
    VoxelToolTerrain(VoxelTerrain *terrain, VoxelMap* map);

    bool is_area_editable(const Rect3i &box) const override;
    VoxelRaycastResult* raycast(Vector3 pos, Vector3 dir, float max_distance, uint32_t collision_mask) override { return nullptr; };

protected:
    uint64_t _get_voxel(Vector3i pos) override;
    float _get_voxel_f(Vector3i pos) override;
    void _set_voxel(Vector3i pos, uint64_t v) override;
    void _set_voxel_f(Vector3i pos, float v) override;
    void _post_edit(const Rect3i &box) override;

    uint64_t _get_voxel(Vector3i pos, int channel) override;
    float _get_voxel_f(Vector3i pos, int channel) override;
    void _set_voxel(Vector3i pos, uint64_t v, int channel) override;
    void _set_voxel_f(Vector3i pos, float v, int channel) override;
private:
    VoxelTerrain *_terrain = nullptr;
    VoxelMap* _map;
};



#endif // VOXELTOOLTERRAIN_H
