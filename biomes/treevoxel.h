#ifndef TREEVOXEL_H
#define TREEVOXEL_H

#include "commonMath/math_funcs.h"
#include "commonMath/cone.h"
#include "commonMath/vector2.h"

struct RoundCone {
    Vector3 m_vOrigin;
    float m_fRadiusDown, m_fRadiusUp;
    float m_fHeight;

    RoundCone(Vector3 vOrigin, float fRadiusDown, float fRadiusUp, float fHeight)
        : m_vOrigin(vOrigin), m_fRadiusDown(fRadiusDown), m_fRadiusUp(fRadiusUp), m_fHeight(fHeight) {};

    Boxi getBox() {
        float fmaxRadius = std::max(m_fRadiusDown, m_fRadiusUp);
        return Boxi(m_vOrigin - Vector3(fmaxRadius, m_fRadiusDown, fmaxRadius), m_vOrigin + Vector3(fmaxRadius, m_fHeight + m_fRadiusUp, fmaxRadius));
    }
};

class TreeVoxel
{
public:
    enum Ease_Type {
        Ease_linear,
        Ease_in_quad,
        Ease_out_quad
    };

    enum class OBJECT_TYPES : uint8_t {
            Shrub = 0, // 灌木
            Arbor, // 乔木
            Lianas, // 藤本
            Herb, // 草本
            Tree1,
            Tree2,
            Tree3,
            Tree4
        };


    TreeVoxel();

    //static void setTree(VoxelToolTerrain* pVoxelTool, const Vector3& vroot, float fradius, float fheight);

    static float capped_cone(const Vector3i &posi, Cone cone, float ra, float rb);
    static float rounded_cone(const Vector3& pos, const Vector3& origin, float fradius_min, float fradius_max, float fheight);
    static float capsule(const Vector3& pos, Vector3 vposa, Vector3 vposb, float fradius);
    static Vector3 bend_linear(const Vector3& pos, Vector3 vposa, Vector3 vposb, Vector3 vpos);
    static float slab();
    static float torus(const Vector3& pos, float fradiusa, float fradiusb);
    static float image(const Vector3i &posi, std::string qsImageName);

    static float ease(float v, Ease_Type eType = Ease_linear);
};

#endif // TREEVOXEL_H
