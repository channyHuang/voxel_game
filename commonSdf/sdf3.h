#ifndef SDF3_H
#define SDF3_H

#include <functional>
#include <vector>
using namespace std;

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

    void Sphere(real fradius, Vector3 vcenter);
    std::function<real(const Vector3 vpos)> Plane(Vector3 vnormal, Vector3 vpoint);
    std::function<real(const Vector3 vpos)> Slab(Vector3 v0, Vector3 v1, float k = 0);
    std::function<real(const Vector3 vpos)> Box(Vector3 vsize, Vector3 vcenter);
    std::function<real(const Vector3 vpos)> RoundBox(Vector3 vsize, real fradius);
    std::function<real(const Vector3 vpos)> WireframeBox(Vector3 vsize, real fthickness);
    std::function<real(const Vector3 vpos)> Torus(real fradius1, real fradius2);
    std::function<real(const Vector3 vpos)> Capsule(Vector3 vposa, Vector3 vposb, real fradius);
    std::function<real(const Vector3 vpos)> Cylinder(real fradius);
    std::function<real(const Vector3 vpos)> CappedCylinder(Vector3 vposa, Vector3 vposb, real fradius);
    std::function<real(const Vector3 vpos)> RoundedCylinder(real fradiusa, real fradiusb, real height);
    std::function<real(const Vector3 vpos)> CappedCone(Vector3 vposa, Vector3 vposb, real fradiusa, real fradiusb);
    std::function<real(const Vector3 vpos)> RoundedCone(real fradius1, real fradius2, real height);
    std::function<real(const Vector3 vpos)> Ellipsoid(real fsize);
    std::function<real(const Vector3 vpos)> Pyramid(real height);
    std::function<real(const Vector3 vpos)> Tetrahedron(real fradius);
    std::function<real(const Vector3 vpos)> Octahedron(real fradius);
    std::function<real(const Vector3 vpos)> Dodecahedron(real fradius);
    std::function<real(const Vector3 vpos)> Icosahedron(real fradius);

    std::function<real(const Vector3 vpos)> TriPrism(Vector2 height);

    // --------------------------------------------------------------
    // operation
    // --------------------------------------------------------------

    std::function<real(const Vector3 vpos)> translate(std::function<real(const Vector3 vpos)> funSdf, Vector3 voffset);
    std::function<real(const Vector3 vpos)> scale(std::function<real(const Vector3 vpos)> funSdf, Vector3 vdir);
    std::function<real(const Vector3 vpos)> rotate(std::function<real(const Vector3 vpos)> funSdf, real fangle, Vector3 vdir);
    std::function<real(const Vector3 vpos)> rotate_to(std::function<real(const Vector3 vpos)> funSdf, Vector3 va, Vector3 vb);
    std::function<real(const Vector3 vpos)> orient(std::function<real(const Vector3 vpos)> funSdf, Vector3 vaxis);
    std::function<real(const Vector3 vpos)> circular_array(std::function<real(const Vector3 vpos)> funSdf, int count);
    std::function<real(const Vector3 vpos)> elongate(std::function<real(const Vector3 vpos)> funSdf, Vector3 vsize);
    std::function<real(const Vector3 vpos)> twist(std::function<real(const Vector3 vpos)> funSdf, int k);
    std::function<real(const Vector3 vpos)> bend(std::function<real(const Vector3 vpos)> funSdf, int k);
    std::function<real(const Vector3 vpos)> bend_linear(std::function<real(const Vector3 vpos)> funSdf, Vector3 vpa, Vector3 vpb, Vector3 v);
    std::function<real(const Vector3 vpos)> bend_radial(std::function<real(const Vector3 vpos)> funSdf);
    std::function<real(const Vector3 vpos)> transition_linear(std::function<real(const Vector3 vpos)> funSdf);
    std::function<real(const Vector3 vpos)> transition_radial(std::function<real(const Vector3 vpos)> funSdf);
    std::function<real(const Vector3 vpos)> wrap_around(std::function<real(const Vector3 vpos)> funSdf);
    std::function<real(const Vector3 vpos)> slice(std::function<real(const Vector3 vpos)> funSdf);

    // --------------------------------------------------------------
    // dn
    // --------------------------------------------------------------

    std::function<real(const Vector3 vpos)> Union(vector<std::function<real(const Vector3 vpos)> > funs);
    std::function<real(const Vector3 vpos)> difference(vector<std::function<real(const Vector3 vpos)> > funs);
    std::function<real(const Vector3 vpos)> intersection(vector<std::function<real(const Vector3 vpos)> > funs, float k);
    std::function<real(const Vector3 vpos)> blend(vector<std::function<real(const Vector3 vpos)> > funs);
    std::function<real(const Vector3 vpos)> negate(vector<std::function<real(const Vector3 vpos)> > funs);
    std::function<real(const Vector3 vpos)> dilate(vector<std::function<real(const Vector3 vpos)> > funs);
    std::function<real(const Vector3 vpos)> erode(std::function<real(const Vector3 vpos)> funSdf, real r);
    std::function<real(const Vector3 vpos)> shell(std::function<real(const Vector3 vpos)> funSdf, real r);
    std::function<real(const Vector3 vpos)> repeat(std::function<real(const Vector3 vpos)> funSdf, real spacing);

private:
    std::function<real(const Vector3 vpos)> funSdf;
};

#endif // SDF3_H
