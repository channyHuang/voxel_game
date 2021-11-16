#include "treevoxel.h"


TreeVoxel::TreeVoxel()
{

}

void TreeVoxel::setTree(VoxelToolTerrain* pVoxelTool, const Vector3& vroot, float fradius, float fheight) {
    // brench
    for (int i = 0; i <= fradius + 1; ++i) {
        pVoxelTool->set_voxel_info(vroot + Vector3(0, i, 0), -1, 1);
    }
    Cone cone(vroot + Vector3(0, fradius - 1 , 0), fradius, fheight);
    Boxi boxi = cone.getBox();

    Vector3 vpos;
    for (vpos.x = boxi.vMin.x; vpos.x < boxi.vMax.x; vpos.x++) {
        for (vpos.y = boxi.vMin.y; vpos.y < boxi.vMax.y; vpos.y++) {
            for (vpos.z = boxi.vMin.z; vpos.z < boxi.vMax.z; vpos.z++) {
                float sdf = capped_cone(vpos, cone, 0, 5);
                if (sdf < 0)
                      pVoxelTool->set_voxel_info(vpos, sdf, (sdf > 0 ? 0 : 1));
            }
        }
    }

}

float TreeVoxel::capped_cone(const Vector3i &posi, Cone cone, float ra, float rb) {
    Vector3i a = cone.getCenter(), b = cone.getTop();
    float rba = rb - ra;
    float baba = (b - a).dot(b - a);
    float papa = (posi - a).dot(posi - a);
    float paba = (posi - a).dot(b - a) / baba;
    float x = std::sqrt(papa - paba * paba * baba);
    float cax = std::max(0.f, x - (paba < 0.5 ? ra : rb));
    float cay = std::abs(paba - 0.5f) - 0.5f;
    float k = rba * rba + baba;
    float f = Math::clamp((rba * (x - ra) + paba * baba) / k, 0.f, 1.f);
    float cbx = x - ra - f * rba;
    float cby = paba - f;
    float s = (cbx < 0 && cay < 0 ? -1 : 1);
    return s * std::sqrt(std::min(cax * cax + cay * cay * baba, cbx * cbx + cby * cby * baba));
}
