#ifndef VOXEL_MESHER_H
#define VOXEL_MESHER_H

#include "../voxels/voxel_buffer.h"
#include "../common_math/vector3.h"

#include <array>

struct Arrays {
    std::vector<Vector3> positions;
    std::vector<Vector3> normals;
    std::vector<uint32_t> indices;
    bool isWater;

    bool empty() {
        return (positions.size() <= 0 || indices.size() <= 0);
    }
};

class VoxelMesher {
public:
	struct Input {
		const VoxelBuffer &voxels;
        int lod;
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
	int _minimum_padding = 0;
	int _maximum_padding = 0;
};

#endif // VOXEL_MESHER_H
