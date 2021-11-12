#include "voxeltool.h"


Vector3 VoxelRaycastResult::_b_get_position() const {
    return position.to_vec3();
}

Vector3 VoxelRaycastResult::_b_get_previous_position() const {
    return previous_position.to_vec3();
}

//----------------------------------------

void VoxelTool::set_value(uint64_t val) {
    _value = val;
}

uint64_t VoxelTool::get_value() const {
    return _value;
}

void VoxelTool::set_eraser_value(uint64_t value) {
    _eraser_value = value;
}

uint64_t VoxelTool::get_eraser_value() const {
    return _eraser_value;
}

void VoxelTool::set_channel(int channel) {
    _channel = channel;
}

int VoxelTool::get_channel() const {
    return _channel;
}

void VoxelTool::set_mode(Mode mode) {
    _mode = mode;
}

VoxelTool::Mode VoxelTool::get_mode() const {
    return _mode;
}

VoxelRaycastResult* VoxelTool::raycast(Vector3 pos, Vector3 dir, float max_distance, uint32_t collision_mask) {
    return new VoxelRaycastResult();
}

uint64_t VoxelTool::get_voxel(Vector3i pos) {
    return _get_voxel(pos);
}

uint64_t VoxelTool::get_voxel(Vector3i pos, int channel) {
    return _get_voxel(pos, channel);
}

float VoxelTool::get_voxel_f(Vector3i pos) {
    return _get_voxel_f(pos);
}

float VoxelTool::get_voxel_f(Vector3i pos, int channel) {
    return _get_voxel_f(pos, channel);
}

void VoxelTool::set_voxel(Vector3i pos, uint64_t v) {
    Rect3i box(pos, Vector3i(1));
    if (!is_area_editable(box)) {
        return;
    }
    _set_voxel(pos, v);
    _post_edit(box);
}

void VoxelTool::set_voxel(Vector3i pos, uint64_t v, int channel) {
    Rect3i box(pos, Vector3i(1));
    if (!is_area_editable(box)) {
        return;
    }
    _set_voxel(pos, v, channel);
    _post_edit(box);
}

void VoxelTool::set_voxel_f(Vector3i pos, float v) {
    Rect3i box(pos, Vector3i(1));
    if (!is_area_editable(box)) {
        return;
    }
    _set_voxel_f(pos, v);
    _post_edit(box);
}

void VoxelTool::set_voxel_f(Vector3i pos, float v, int channel) {
    Rect3i box(pos, Vector3i(1));
    if (!is_area_editable(box)) {
        return;
    }
    _set_voxel_f(pos, v, channel);
    _post_edit(box);
}

void VoxelTool::do_point(Vector3i pos) {
    Rect3i box(pos, Vector3i(1));
    if (!is_area_editable(box)) {
        return;
    }
    if (_channel == VoxelBuffer::CHANNEL_SDF) {
        _set_voxel_f(pos, _mode == MODE_REMOVE ? 1.0 : -1.0);
    } else {
        _set_voxel(pos, _mode == MODE_REMOVE ? _eraser_value : _value);
    }
    _post_edit(box);
}

void VoxelTool::do_line(Vector3i begin, Vector3i end) {
}

void VoxelTool::do_circle(Vector3i pos, int radius, Vector3i direction) {
}

uint64_t VoxelTool::_get_voxel(Vector3i pos) {
    return 0;
}

uint64_t VoxelTool::_get_voxel(Vector3i pos, int channel) {
    return 0;
}

float VoxelTool::_get_voxel_f(Vector3i pos) {
    return 0;
}

float VoxelTool::_get_voxel_f(Vector3i pos, int channel) {
    return 0;
}

void VoxelTool::_set_voxel(Vector3i pos, uint64_t v) {
}

void VoxelTool::_set_voxel(Vector3i pos, uint64_t v, int channel) {
}


void VoxelTool::_set_voxel_f(Vector3i pos, float v) {
}

void VoxelTool::_set_voxel_f(Vector3i pos, float v, int channel) {
}

namespace {
inline float sdf_blend(float src_value, float dst_value, VoxelTool::Mode mode) {
    float res;
    switch (mode) {
        case VoxelTool::MODE_ADD:
            // Union
            res = std::min(src_value, dst_value);
            break;

        case VoxelTool::MODE_REMOVE:
            // Relative complement (or difference)
            res = std::max(1.f - src_value, dst_value);
            break;

        case VoxelTool::MODE_SET:
            res = src_value;
            break;

        default:
            res = 0;
            break;
    }
    return res;
}
} // namespace


bool VoxelTool::is_area_editable(const Rect3i &box) const {
    return false;
}

void VoxelTool::_post_edit(const Rect3i &box) {
}


