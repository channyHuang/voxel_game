#pragma once

#include "vector3.h"

class Box
{
public:
    Box() { vMin = Vector3(0), vMax = Vector3(0); };
    Box(const Vector3 &_vMin, const Vector3 &_vMax) : vMin(_vMin), vMax(_vMax) {}
    ~Box() {};

public:
    Vector3 vMin, vMax;
};
