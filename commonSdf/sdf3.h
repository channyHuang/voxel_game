#ifndef SDF3_H
#define SDF3_H

#include <functional>

#include "vector2.h"
#include "vector3.h"
#include "math_funcs.h"
#include "typeheader.h"

class Sdf3
{
public:
    Sdf3();

    void reset();

    // --------------------------------------------------------------
    // calculation
    // --------------------------------------------------------------

    void Sphere(real fradius,  Vector3 vcenter);
    void Plane(Vector3 vnormal, Vector3 vpoint);
    void Slab(Vector3 v0, Vector3 v1);
    void Box(Vector3 vsize, Vector3 vcenter);
    void RoundBox(Vector3 vsize, real fradius);
    void WireframeBox(Vector3 vsize, real fthickness);
    void Torus(real fradius1, real fradius2);
    void Capsule(Vector3 vposa, Vector3 vposb, real fradius);
    void Cylinder(real fradius);
    void CappedCylinder(Vector3 vposa, Vector3 vposb, real fradius);
    void RoundedCylinder(real fradiusa, real fradiusb, real height);
    void CappedCone(Vector3 vposa, Vector3 vposb, real fradiusa, real fradiusb);
    void RoundedCone(real fradius1, real fradius2, real height);
    void Ellipsoid(real fsize);
    void Pyramid(real height);
    void Tetrahedron(real fradius);
    void Octahedron(real fradius);
    void Dodecahedron(real fradius);
    void Icosahedron(real fradius);

    void TriPrism(Vector2 height);
    // --------------------------------------------------------------
    // operation
    // --------------------------------------------------------------

    void translate(std::function<real(const Vector3 vpos)> fun, Vector3 voffset);
    void scale(std::function<real(const Vector3 vpos)> fun, Vector3 vdir);
    void rotate(std::function<real(const Vector3 vpos)> fun, real fangle, Vector3 vdir);

private:
    std::function<real(const Vector3 vpos)> funSdf;
};

#endif // SDF3_H
