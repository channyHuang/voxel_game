#include "treevoxel.h"

#include "common_math/cone.h"

TreeVoxel::TreeVoxel()
{

}

void TreeVoxel::setTree(VoxelToolTerrain* pVoxelTool, const Vector3& vroot, float fradius, float fheight) {
    // brench
    for (int i = 1; i <= fradius; ++i) {
        pVoxelTool->set_voxel_info(vroot + Vector3(0, i, 0), -1, 1);
    }
    Cone cone(vroot + Vector3(0, fradius , 0), fradius, fheight);
    Boxi boxi = cone.getBox();

    Vector3 vpos;
    for (vpos.x = boxi.vMin.x; vpos.x < boxi.vMax.x; vpos.x++) {
        for (vpos.y = boxi.vMin.y; vpos.y < boxi.vMax.y; vpos.y++) {
            for (vpos.z = boxi.vMin.z; vpos.z < boxi.vMax.z; vpos.z++) {
                float sdf = cone.getSdf(vpos);
                if (sdf <= 0)
                    pVoxelTool->set_voxel_info(vpos, sdf, (sdf > 0 ? 0 : 1));
            }
        }
    }

}
