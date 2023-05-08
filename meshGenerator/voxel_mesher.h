#ifndef VOXEL_MESHER_H
#define VOXEL_MESHER_H

#include <array>

#include "commonMath/vector3.h"

#include "terrains/terrainCommonStruct.h"

class VoxelMesher {
public:
    virtual void build(MeshOutput &output, const MeshInput &voxels);

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
