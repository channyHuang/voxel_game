#include "sdf3.h"

Sdf3::Sdf3()
{
    funSdf = nullptr;
}

void Sdf3::reset() {
    funSdf = nullptr;
}

// --------------------------------------------------------------
// calculation
// --------------------------------------------------------------

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

void Sdf3::Slab(Vector3 v0, Vector3 v1) {
    // not impliment yeah
}

void Sdf3::Box(Vector3 vsize, Vector3 vcenter) {
    funSdf = [&](const Vector3 vpos) -> real {
        Vector3 vdist = (vpos - vcenter) - vsize * 0.5f;
        return (vdist.getMax(0)).len() + std::min(vdist.getMaxNum(), 0.f);
    };
}

void Sdf3::RoundBox(Vector3 vsize, real fradius) {
    funSdf = [&](const Vector3 vpos) -> real {
        Vector3 vdist = vpos.getAbs() - vsize * 0.5f + fradius;
        return (vdist.getMax(0)).len() + std::fmin(vdist.getMaxNum(), 0) - fradius;
    };
}

void Sdf3::WireframeBox(Vector3 vsize, real fthickness) {
    std::function<real(const Vector3 v)> funG = [&](Vector3 v) -> real {
        return Vector3(std::fmax(v.x, 0), std::fmax(v.y, 0), std::fmax(v.z, 0)).len() + std::fmin(std::fmax(v.x, std::fmax(v.y, v.z)), 0);
    };

    funSdf = [&](const Vector3 vpos) -> real {
        Vector3 vp = vpos.getAbs() - vsize / 2 - fthickness / 2;
        Vector3 vq = (vp + fthickness / 2).getAbs() - fthickness / 2;
        return std::fmin(std::fmin(funG(Vector3(vp.x, vq.y, vq.z)), funG(Vector3(vq.x, vp.y, vq.z))), funG(Vector3(vq.x, vq.y, vp.z)));
    };
}

void Sdf3::Torus(real fradius1, real fradius2) {
    funSdf = [&](const Vector3 vpos) -> real {
        real a = Vector2(vpos.x, vpos.y).len() - fradius1;
        return Vector2(a, vpos.z).len() - fradius2;
    };
}

void Sdf3::Capsule(Vector3 vposa, Vector3 vposb, real fradius) {
    funSdf = [&](const Vector3 vpos) -> real {
        Vector3 vpa = vpos - vposa;
        Vector3 vba = vposb - vposa;
        real h = Math::Clamp(vpa.dot(vba) / vba.dot(vba), 0.f, 1.f);
        return (vpa - vba * h).len() - fradius;
    };
}

void Sdf3::Cylinder(real fradius) {
    funSdf = [&](const Vector3 vpos) -> real {
        return Vector2(vpos.x, vpos.y).len() - fradius;
    };
}

void Sdf3::CappedCylinder(Vector3 vposa, Vector3 vposb, real fradius) {
    funSdf = [&](const Vector3 vpos) -> real {
        Vector3 vba = vposb - vposa;
        real baba = vba.dot(vba);
        real paba = vposa.dot(vba);
        real x = (vposa * baba - vba * paba).len() - fradius * baba;
        real y = std::fabs(paba - baba * 0.5f) - baba * 0.5f;
        real d = (std::fmax(x, y) < 0 ? -std::fmin(x * x, y * y * baba)
                                      : (x > 0 ? x * x : 0) + (y > 0 ? y * y * baba : 0));
        return std::signbit(d) * std::sqrt(std::fabs(d)) / baba;
    };
}

void Sdf3::RoundedCylinder(real fradiusa, real fradiusb, real height) {
    funSdf = [&](const Vector3 vpos) -> real {
        Vector2 vd = Vector2(Vector2(vpos.x, vpos.y).len() - fradiusa + fradiusb, std::fabs(vpos.z) - height / 2.f + fradiusb);
        return std::fmin(std::fmax(vd.x, vd.y), 0) + (Vector2(std::fmax(vd.x, 0), std::fmax(vd.y, 0)).len() - fradiusb);
    };
}

void Sdf3::CappedCone(Vector3 vposa, Vector3 vposb, real fradiusa, real fradiusb) {
    funSdf = [&](const Vector3 vpos) -> real {
        real frba = fradiusa - fradiusb;
        real fbaba = (vposb - vposa).dot(vposb - vposa);
        real fpapa = (vpos - vposa).dot(vpos - vposa);
        real fpaba = (vpos - vposa).dot(vposb - vposa) / fbaba;
        real fx = std::sqrt(fpapa - fpaba * fpaba * fbaba);
        real fcax = std::fmax(0, fx - (fpaba < 0.5f ? fradiusa : fradiusb));
        real fk = frba * frba + fbaba;
        real ff = Math::Clamp((frba * (fx - fradiusa) + fpaba * fbaba) / fk, 0.f, 1.f);
        real fcay = std::fabs(fpaba - 0.5) - 0.5f;
        real fcbx = fx - fradiusa - ff * frba;
        real fcby = fpaba - ff;
        real fs = ((fcbx < 0 && fcay < 0) ? -1 : 1);
        return fs * std::sqrt(std::fmin(fcax * fcax + fcay * fcay * fbaba, fcbx * fcbx + fcby * fcby * fbaba));
     };
}

void Sdf3::RoundedCone(real fradius1, real fradius2, real height) {
    funSdf = [&](const Vector3 vpos) -> real {
        Vector2 vq = Vector2(Vector2(vpos.x, vpos.y).len(), vpos.z);
        real fb = (fradius1 - fradius2) / height;
        real fa = std::sqrt(1 - fb * fb);
        real fk = vq.dot(Vector2(-fb, fa));
        real c1 = vq.len() - fradius1;
        real c2 = (vq - Vector2(0, height)).len() - fradius2;
        real c3 = vq.dot(Vector2(fa, fb)) - fradius1;
        return (fk < 0 ? c1 : (fk > fa * height ? c2 : c3));
    };
}

void Sdf3::Ellipsoid(real fsize) {
    funSdf = [&](const Vector3 vpos) -> real {
        real k0 = (vpos / fsize).len();
        real k1 = (vpos / (fsize * fsize)).len();
        return k0 * (k0 - 1) / k1;
    };
}

void Sdf3::Pyramid(real height) {
    funSdf = [&](const Vector3 vpos) -> real {
        Vector2 va = Vector2(vpos.x, vpos.y).getAbs() - 0.5f;
        int w = (va.y > va.x);
        va[w] = Vector2(va[1], va[0])[w];
        Vector3 vp = Vector3(va.x, vpos.z, va.y);
        real m2 = height * height + 0.25;
        Vector3 vq = Vector3(vp.z, height * vp.y - 0.5f * vp.x, height * vp.x + 0.5f * vp.y);
        real fs = std::fmax(-vq.x, 0);
        real ft = Math::Clamp((vq.y - 0.5f * vp.z) / (m2 + 0.25f), 0.f, 1.f);
        real fa = m2 * std::pow((vq.x + fs), 2) + vq.y * vq.y;
        real fb = m2 * std::pow(vq.x + 0.5f * ft, 2) + std::pow(vq.y - m2 * ft, 2);
        real d2 = (std::fmin(vq.y, -vq.x * m2 - vq.y * 0.5f) > 0 ? 0 : std::fmin(fa, fb));
        return std::sqrt((d2 + vq.z * vq.z) / m2) * std::signbit(std::fmax(vq.z, -vp.y));
    };
}

void Sdf3::Tetrahedron(real fradius) {
    funSdf = [&](const Vector3 vpos) -> real {
        return (std::max(std::fabs(vpos.x + vpos.y) - vpos.z, std::fabs(vpos.x - vpos.y) + vpos.z) - 1) / std::sqrt(3);
    };
}

void Sdf3::Octahedron(real fradius) {
    funSdf = [&](const Vector3 vpos) -> real {
        return (vpos.getAbs().getSumNum() - fradius) * std::tan(Math::PI / 60);
    };
}

void Sdf3::Dodecahedron(real fradius) {
    Vector3 v = Vector3((std::sqrt(5) + 1) / 2.f, 1.f, 0.f);
    funSdf = [&](const Vector3 vpos) -> real {
        Vector3 vp = (vpos / fradius).getAbs();
        real a = vp.dot(v);
        real b = vp.dot(Vector3(v.z, v.x, v.y));
        real c = vp.dot(Vector3(v.y, v.z, v.x));
        real q = (std::fmax(std::fmax(a, b), c) - v.x) * fradius;
        return q;
    };
}

void Sdf3::Icosahedron(real fradius) {
    Vector3 v = Vector3((std::sqrt(5) + 3) / 2.f, 1.f, 0.f);
    real w = std::sqrt(3) / 3;
    funSdf = [&](const Vector3 vpos) -> real {
        Vector3 vp = (vpos / fradius).getAbs();
        real a = vp.dot(v);
        real b = vp.dot(Vector3(v.z, v.x, v.y));
        real c = vp.dot(Vector3(v.y, v.z, v.x));
        real d = vp.dot(Vector3(w));
        return std::fmax(std::fmax(std::fmax(a, b), c) - v.x, d) * fradius;
    };
}

void Sdf3::TriPrism(Vector2 height) {
    funSdf = [&](const Vector3 vpos) -> real {
        Vector3 vq = vpos.getAbs();
        return std::fmax(vq.z - height.y, std::fmax(vq.x * 0.866025 + vpos.y * 0.5, -vpos.y) - height.x * 0.5f);
    };
}

// --------------------------------------------------------------
// operation
// --------------------------------------------------------------

void Sdf3::translate(std::function<real(const Vector3 vpos)> fun, Vector3 voffset) {
    funSdf = [&](const Vector3 vpos) -> real {
        return fun(vpos - voffset);
    };
}

void Sdf3::scale(std::function<real(const Vector3 vpos)> fun, Vector3 vdir) {
    // no idea
}

void Sdf3::rotate(std::function<real(const Vector3 vpos)> fun, real fangle, Vector3 vdir) {
    funSdf = [&](const Vector3 vpos) -> real {
        //return fun(vpos.dot());
    };
}


