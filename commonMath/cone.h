#ifndef CONE_H
#define CONE_H

#include "vector3.h"
#include "boxi.h"

class Cone
{
public:
    Cone(Vector3 vcenter, float fRadius, float fHeight)
        : m_vCenter(vcenter), m_fRadius(fRadius), m_fHeight(fHeight) {};

    float getSdf(Vector3 vPos);
    Boxi getBox();
    Vector3 getCenter() { return m_vCenter; }
    Vector3 getTop() { return m_vCenter + Vector3(0, m_fHeight, 0); }

    static Vector3i vector3FloorOrCeil(const Vector3 &pos, bool bfloor = true) {
        if (bfloor) {
            Vector3 &&pos2floor = pos.getFloor();
            Vector3i &&vposi = Vector3i(static_cast<int>(pos2floor.x), static_cast<int>(pos2floor.y), static_cast<int>(pos2floor.z));
            return vposi;
        }
        //bceil
        Vector3 &&pos2ceil = Vector3(pos);
        pos2ceil.ceil();
        Vector3i &&vposi = Vector3i(static_cast<int>(pos2ceil.x), static_cast<int>(pos2ceil.y), static_cast<int>(pos2ceil.z));
        return vposi;
    }
private:
    Vector3 m_vCenter;
    float m_fRadius;
    float m_fHeight;
};

#endif // CONE_H
