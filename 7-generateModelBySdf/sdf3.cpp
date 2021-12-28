#include "sdf3.h"

Sdf3::Sdf3()
{
    funSdf = nullptr;
}

void Sdf3::reset() {
    funSdf = nullptr;
}

void Sdf3::Sphere(real fradius, const Vector3 vcenter) {
    funSdf = [&](const Vector3 vpos) -> real {
        return (vpos - vcenter).len() - fradius;
    };
}

void Sdf3::Plane(const Vector3 vnormal, Vector3 vpoint) {
    Vector3 vnorm = vnormal.getNormalize();
    funSdf = [&](const Vector3 vpos) -> real {
        return (vpoint - vpos).dot(vnorm);
    };
}
