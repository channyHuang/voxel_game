#ifndef TREEVOXEL_H
#define TREEVOXEL_H

#include "terrains/voxeltoolterrain.h"
#include "commonMath/math_funcs.h"
#include "commonMath/cone.h"


class TreeVoxel
{
public:
    TreeVoxel();

    static void setTree(VoxelToolTerrain* pVoxelTool, const Vector3& vroot, float fradius, float fheight);

    static float capped_cone(const Vector3i &posi, Cone cone, float ra, float rb);
};

#endif // TREEVOXEL_H
