#pragma once

#include "vector3i.h"

class Boxi
{
public:
    Boxi() {};
    Boxi(const Vector3i &_vMin, const Vector3i &_vMax) : vMin(_vMin), vMax(_vMax) {}
    ~Boxi() {};

    inline bool operator == (const Boxi &v) const {
        return (vMin == v.vMin && vMax == v.vMax);
    }

    bool contains(const Boxi &box) {
        return (vMin.x < box.vMin.x && vMin.y < box.vMin.y
                && vMax.x > box.vMax.x && vMax.y > box.vMax.y);
    };

    bool contains(const Vector3i &pos) {
        return (vMin.x < pos.x && vMin.y < pos.y && vMin.z < pos.z && vMax.x > pos.x && vMax.y > pos.y && vMax.z > pos.z);
    }

    static Boxi from_center_extents(const Vector3i &_vCenter, const Vector3i &_vExpend) {
        return Boxi(_vCenter - _vExpend, _vCenter + _vExpend);
    }

    bool intersects(const Boxi& b) const {
        if (vMin.x > b.vMax.x || vMin.y > b.vMax.y || vMin.z > b.vMax.z) return false;
        if (vMax.x < b.vMin.x || vMax.y < b.vMin.y || vMax.z < b.vMax.z) return false;
        return true;
    }

    template<typename A>
    void difference(const Boxi &b, A action) {
        if (!intersects(b)) {
            action(*this);
            return;
        }
        Boxi target = *this;

        if (target.vMin.x < b.vMin.x) {
            action(Boxi(target.vMin, Vector3i(b.vMin.x, target.vMax.y, target.vMax.z)));
            target.vMin.x = b.vMin.x;
        }
        if (target.vMin.y < b.vMin.y) {
            action(Boxi(target.vMin, Vector3i(target.vMax.x, b.vMin.y, target.vMax.z)));
            target.vMin.y = b.vMin.y;
        }
        if (target.vMin.z < b.vMin.z) {
            action(Boxi(target.vMin, Vector3i(target.vMax.x, target.vMax.y, b.vMin.z)));
            target.vMin.z = b.vMin.z;
        }
        if (target.vMax.x > b.vMax.x) {
            action(Boxi(Vector3i(b.vMax.x, target.vMin.y, target.vMin.z), target.vMax));
            target.vMax.x = b.vMax.x;
        }
        if (target.vMax.y > b.vMax.y) {
            action(Boxi(Vector3i(target.vMin.x, b.vMax.y, target.vMin.z), target.vMax));
            target.vMax.y = b.vMax.y;
        }
        if (target.vMax.z > b.vMax.z) {
            action(Boxi(Vector3i(target.vMin.x, target.vMin.y, b.vMax.z), target.vMax));
        }
    }
    template<typename A>
    void for_each_cell(A a) const {
        Vector3i posi;
        for (posi.x = vMin.x; posi.x < vMax.x; ++posi.x) {
            for (posi.y = vMin.y; posi.y < vMax.y; ++posi.y) {
                for (posi.z = vMin.z; posi.z < vMax.z; ++posi.z) {
                    a(posi);
                }
            }
        }
    }

public:
    Vector3i vMin, vMax;
};


