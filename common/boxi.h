#ifndef BOXI_H
#define BOXI_H

#include "vector3i.h"

class Boxi
{
public:
    Boxi();
    Boxi(const Vector3i &_vMin, const Vector3i &_vMax) : vMin(_vMin), vMax(_vMax) {}
    ~Boxi();

public:
    Vector3i vMin, vMax;
};

#endif // BOXI_H
