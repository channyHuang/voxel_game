#include "treevoxel.h"
#include "common_enum.h"

TreeVoxel::TreeVoxel()
{

}

void TreeVoxel::setTree(VoxelToolTerrain* pVoxelTool, const Vector3& vroot, float fradius, float fheight) {
        float sdf = -std::max(0.5f, rand() / (float)RAND_MAX), fold_sdf;
        Vector3 vpos;
        Boxi box;

        OBJECT_TYPES eobject_type = (OBJECT_TYPES)(OBJECT_TYPES::Herb);
        switch (eobject_type) {
        case OBJECT_TYPES::Shrub:
        {
            // brench
            for (int i = 0; i <= 3; ++i) {
                pVoxelTool->set_voxel_info(vroot + Vector3(0, i, 0), sdf, MaterialType::MUD);
            }
            Vector3 vorigin = vroot + Vector3(0, 2, 0) + Vector3(0, 2, 0);
            RoundCone cone(vorigin, 2, 0, 2);
            box = cone.getBox();
            for (vpos.x = box.vMin.x; vpos.x <= box.vMax.x; vpos.x++) {
                for (vpos.y = box.vMin.y; vpos.y <= box.vMax.y; vpos.y++) {
                    for (vpos.z = box.vMin.z; vpos.z <= box.vMax.z; vpos.z++) {
                        sdf = rounded_cone(vpos, vorigin, cone.m_fRadiusDown, cone.m_fRadiusUp, cone.m_fHeight);
                        fold_sdf = pVoxelTool->get_voxel_f(vpos);
                        if (sdf <= 0 && sdf < fold_sdf) {
                            pVoxelTool->set_voxel_info(vpos, sdf, MaterialType::GRASS);
                        }
                    }
                }
            }
            break;
        }
        case OBJECT_TYPES::Arbor:
        {
            // brench
            for (int i = 0; i <= 3; ++i) {
                pVoxelTool->set_voxel_info(vroot + Vector3(0, i, 0), sdf, MaterialType::MUD);
            }
            Vector3 vorigin = vroot + Vector3(0, 3, 0) + Vector3(0, 3, 0);
            RoundCone cone(vorigin, 3, 0, 6);
            box = cone.getBox();
            for (vpos.x = box.vMin.x; vpos.x <= box.vMax.x; vpos.x++) {
                for (vpos.y = box.vMin.y; vpos.y <= box.vMax.y; vpos.y++) {
                    for (vpos.z = box.vMin.z; vpos.z <= box.vMax.z; vpos.z++) {
                        sdf = rounded_cone(vpos, vorigin, cone.m_fRadiusDown, cone.m_fRadiusUp, cone.m_fHeight);
                        fold_sdf = pVoxelTool->get_voxel_f(vpos);
                        if (sdf <= 0 && sdf < fold_sdf) {
                            pVoxelTool->set_voxel_info(vpos, sdf, MaterialType::GRASS);
                        }
                    }
                }
            }
            break;
        }
        case OBJECT_TYPES::Herb:
        {

            box = Boxi(vroot - Vector3(100), vroot + Vector3(100));
            for (vpos.x = box.vMin.x; vpos.x <= box.vMax.x; vpos.x++) {
                for (vpos.y = box.vMin.y; vpos.y <= box.vMax.y; vpos.y++) {
                    for (vpos.z = box.vMin.z; vpos.z <= box.vMax.z; vpos.z++) {
                        sdf = torus(vpos - vroot, 80, 10);
                        //sdf = capsule(bend_linear(vpos - vroot, Vector3::NEG_UNIT_Y, Vector3::UNIT_Y, Vector3::UNIT_X), Vector3::NEG_UNIT_Y * 2.f, Vector3::UNIT_Y * 2.f, 3);
                        pVoxelTool->set_voxel_info(vpos, sdf, MaterialType::GRASS);
                    }
                }
            }
        }
            break;
        case OBJECT_TYPES::Tree4:
        {
            // brench
            for (int i = 0; i <= fradius + 1; ++i) {
                pVoxelTool->set_voxel_info(vroot + Vector3(0, i, 0), sdf, MaterialType::GRASS);
            }
            Cone cone(vroot + Vector3(0, fradius - 1, 0), fradius, fheight);
            box = cone.getBox();

            for (vpos.x = box.vMin.x; vpos.x <= box.vMax.x; vpos.x++) {
                for (vpos.y = box.vMin.y; vpos.y <= box.vMax.y; vpos.y++) {
                    for (vpos.z = box.vMin.z; vpos.z <= box.vMax.z; vpos.z++) {
                        sdf = capped_cone(vpos, cone, 0, 5);
                        if (sdf <= 0) {
                            pVoxelTool->set_voxel_info(vpos, sdf, MaterialType::GRASS);
                        }
                    }
                }
            }
        }
            break;
        default:
            break;
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

float TreeVoxel::slab() {
    return 0;
}

float TreeVoxel::image(const Vector3i &posi, QString qsImageName) {
    return 0;
}

float TreeVoxel::rounded_cone(const Vector3& pos, const Vector3& origin, float fradius_min, float fradius_max, float fheight) {
    Vector3 vpos = pos - origin;
    Vector2 q = Vector2(Vector2(vpos.x, vpos.z).len(), vpos.y);
    float b = (fradius_min - fradius_max) / fheight;
    float a = std::sqrt(1 - b * b);
    float k = q.dot(Vector2(-b, a));
    float c1 = q.len() - fradius_min;
    float c2 = (q - Vector2(0, fheight)).len() - fradius_max;
    float c3 = q.dot(Vector2(a, b)) - fradius_min;
    return (k < 0 ? c1 : (k > a * fheight ? c2 : c3));
}

Vector3 TreeVoxel::bend_linear(const Vector3& pos, Vector3 vposa, Vector3 vposb, Vector3 vpos) {
    vpos = vpos * (-1.f);
    Vector3 ab = vposb - vposa;
    float t = Math::clamp((pos - vposa).dot(ab) / ab.dot(ab), 0.f, 1.f);
    return (pos + t * vpos);
}

float TreeVoxel::capsule(const Vector3& pos, Vector3 vposa, Vector3 vposb, float fradius) {
    Vector3 pa = pos - vposa;
    Vector3 ba = vposb - vposa;
    float h = Math::clamp(pa.dot(ba) / ba.dot(ba), 0.f, 1.f);
    return (pa - ba * h).len() - fradius;
}


float TreeVoxel::torus(const Vector3& pos, float fradiusa, float fradiusb) {
    Vector2 xy = Vector2(pos.x, pos.y);
    float z = pos.z;
    float a = xy.len() - fradiusa;
    float b = Vector2(a, z).len() - fradiusb;
    return b;
}

float TreeVoxel::ease(float v, Ease_Type eType) {
    switch(eType) {
    case Ease_in_quad:
        return v * v;
    case Ease_out_quad:
        return -v * (v - 2);
    default:
        break;
    }
    return v;
}
