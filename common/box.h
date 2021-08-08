#ifndef BOX_H
#define BOX_H

#include "vector3.h"

class Box
{
public:
    Box() {};
    Box(const Vector3 &_vMin, const Vector3 &_vMax) : vMin(_vMin), vMax(_vMax) {}
    ~Box() {};

public:
    Vector3 vMin, vMax;
};

#endif // BOX_H
