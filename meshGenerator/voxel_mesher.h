#ifndef VOXEL_MESHER_H
#define VOXEL_MESHER_H

#include "../voxels/voxel_buffer.h"
#include "vector3.h"

#include <array>

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

class VoxelMesher {
public:
	struct Input {
		const VoxelBuffer &voxels;
        int lod;
        Vector3i position;
	};

	struct Output {
        std::vector<Arrays> surfaces;
	};

	virtual void build(Output &output, const Input &voxels);

	int get_minimum_padding() const;
	int get_maximum_padding() const;

	virtual VoxelMesher *clone();

protected:
	static void _bind_methods();

	void set_padding(int minimum, int maximum);

private:
    int _minimum_padding = 2;
    int _maximum_padding = 2;
};

#endif // VOXEL_MESHER_H
