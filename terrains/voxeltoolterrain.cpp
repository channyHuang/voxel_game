#include "voxeltoolterrain.h"

#include "terrainCommonStruct.h"

VoxelToolTerrain::VoxelToolTerrain() {
}

VoxelToolTerrain::VoxelToolTerrain(VoxelTerrain *terrain, VoxelMap* map) {
    _terrain = terrain;
    _map = map;
    // Don't destroy the terrain while a voxel tool still references it
}

bool VoxelToolTerrain::is_area_editable(const Rect3i &box) const {
    return true;//_map->is_area_fully_loaded(box.padded(1));
}


uint64_t VoxelToolTerrain::_get_voxel(Vector3i pos) {
    return _map->get_voxel(pos, _channel);
}

float VoxelToolTerrain::_get_voxel_f(Vector3i pos) {
    return _map->get_voxel_f(pos, _channel);
}

uint64_t VoxelToolTerrain::_get_voxel(Vector3i pos, int channel) {
    return _map->get_voxel(pos, channel);
}

float VoxelToolTerrain::_get_voxel_f(Vector3i pos, int channel) {
    return _map->get_voxel_f(pos, channel);
}

void VoxelToolTerrain::_set_voxel(Vector3i pos, uint64_t v) {
    _map->set_voxel(v, pos, _channel);
    _post_edit(Rect3i(pos, pos + 1));
}

void VoxelToolTerrain::_set_voxel(Vector3i pos, uint64_t v, int channel) {
    _map->set_voxel(v, pos, channel);
    _post_edit(Rect3i(pos, pos + 1));
}

void VoxelToolTerrain::_set_voxel_f(Vector3i pos, float v) {
    _map->set_voxel_f(v, pos, _channel);
    _post_edit(Rect3i(pos, pos + 1));
}

void VoxelToolTerrain::_set_voxel_f(Vector3i pos, float v, int channel) {
    _map->set_voxel_f(v, pos, channel);
    _post_edit(Rect3i(pos, pos + 1));
}

void VoxelToolTerrain::_post_edit(const Rect3i &box) {
    //_terrain->make_area_dirty(box);
}



