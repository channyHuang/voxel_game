#ifndef TREEVOXEL_H
#define TREEVOXEL_H

#include "voxeltoolterrain.h"

class TreeVoxel
{
public:
    TreeVoxel();

    static void setTree(VoxelToolTerrain* pVoxelTool, const Vector3& vroot, float fradius, float fheight);


};

#endif // TREEVOXEL_H
