#include "cone.h"
//[-1, 1]
float Cone::getSdf(Vector3 vPos) {
    if (vPos.y <= m_vCenter.y - 1 || vPos.y >= m_vCenter.y + m_fHeight + 1) return 1;
    float fsideLen = std::sqrt(m_fRadius * m_fRadius + m_fHeight * m_fHeight);
    if (vPos.y <= m_vCenter.y) {
        return m_vCenter.y - vPos.y;
    } else if (vPos.y >= m_vCenter.y + m_fHeight) {
        return (m_vCenter + Vector3(0, m_fHeight, 0)).len();
    }
    float fradiusSqr = (vPos.x - m_vCenter.x) * (vPos.x - m_vCenter.x) + (vPos.z - m_vCenter.z) * (vPos.z - m_vCenter.z);
    float fradiusExp = (m_fHeight - vPos.y + m_vCenter.y) * m_fRadius / m_fHeight;

    if (fradiusSqr <= fradiusExp) {
        return -std::min(vPos.y - m_vCenter.y, (fradiusExp - fradiusSqr) * m_fRadius / fsideLen);
    } else {
        return 1.f;
    }
}

Boxi Cone::getBox() {
    Boxi boxi;
    boxi.vMin = vector3FloorOrCeil(m_vCenter - Vector3(m_fRadius, 0, m_fRadius), true);
    boxi.vMax = vector3FloorOrCeil(m_vCenter + Vector3(m_fRadius, m_fHeight, m_fRadius), false);
    return boxi;
}
