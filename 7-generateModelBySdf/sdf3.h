#ifndef SDF3_H
#define SDF3_H

#include "vector3.h"
#include <functional>

#ifndef real
typedef float real;
#endif

class Sdf3
{
public:
    Sdf3();

    void reset();

    void Sphere(real fradius,  Vector3 vcenter);
    void Plane(Vector3 vnormal, Vector3 vpoint);

private:
    std::function<real(const Vector3 vpos)> funSdf;
};

#endif // SDF3_H
