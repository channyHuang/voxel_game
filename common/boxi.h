#pragma once

#include "vector3i.h"

class Boxi
{
public:
    Boxi() {};
    Boxi(const Vector3i &_vMin, const Vector3i &_vMax) : vMin(_vMin), vMax(_vMax) {}
    ~Boxi() {};

    bool contains(const Boxi &box) {
        return (vMin.x < box.vMin.x && vMin.y < box.vMin.y
                && vMax.x > box.vMax.x && vMax.y > box.vMax.y);
    };

public:
    Vector3i vMin, vMax;

};
