#pragma once

#include "vector3.h"

class Triangle {
public:
    Triangle() : bHasNormal(false) {}
    Triangle(const Vector3 &_v0, const Vector3 &_v1, const Vector3 &_v2)
        : v0(_v0), v1(_v1), v2(_v2) {}
    Triangle(const Vector3 &_v0, const Vector3 &_v1, const Vector3 &_v2, const Vector3 &_normal)
        : v0(_v0), v1(_v1), v2(_v2), vNormal(_normal), bHasNormal(true) {}
    ~Triangle() {};

    void calcNormal() {
        if (bHasNormal) return;
        vNormal = (v1 - v0).cross(v2 - v0);
        bHasNormal = true;
    }

    Vector3 getNormal() {
        if (!bHasNormal) {
            calcNormal();
        }
        return vNormal;
    }

    void set(const Vector3 &_v0, const Vector3 &_v1, const Vector3 &_v2) {
        bHasNormal = false;
        v0 = _v0;
        v1 = _v1;
        v2 = _v2;
    }

public:
    Vector3 v0, v1, v2, vNormal;
    bool bHasNormal;
};
