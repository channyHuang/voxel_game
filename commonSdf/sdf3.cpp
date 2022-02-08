#include "sdf3.h"

Sdf3::Sdf3()
{
    fun = nullptr;
}

void Sdf3::reset() {
    fun = nullptr;
}

// --------------------------------------------------------------
// calculation
// --------------------------------------------------------------

std::function<real(const Vector3& vpos)> Sdf3::Sphere(real fradius, const Vector3& vcenter) {
    fun = [fradius, vcenter](const Vector3& vpos) -> real {
        return (vpos - vcenter).len() - fradius;
    };
    return fun;
}

std::function<real(const Vector3& vpos)> Sdf3::Plane(const Vector3 vnormal, Vector3 vpoint) {
    Vector3 vnorm = vnormal.getNormalize();
    fun = [&](const Vector3& vpos) -> real {
        return (vpos - vpoint).dot(vnorm);
    };
    return fun;
}

std::function<real(const Vector3& vpos)> Sdf3::Slab(Vector3 v0, Vector3 v1, float  k) {
    vector<std::function<real(const Vector3& vpos)> > funs;
    funs.push_back(Plane(Vector3(1, 0, 0), Vector3(v0.x, 0, 0)));
    funs.push_back(Plane(Vector3(-1, 0, 0), Vector3(v1.x, 0, 0)));
    funs.push_back(Plane(Vector3(0, 1, 0), Vector3(0, v0.y, 0)));
    funs.push_back(Plane(Vector3(0, -1, 0), Vector3(0, v1.y, 0)));
    funs.push_back(Plane(Vector3(0, 0, 1), Vector3(0, 0, v0.z)));
    funs.push_back(Plane(Vector3(0, 0, -1), Vector3(0, 0, v1.z)));
    return intersection(funs, k);
}

std::function<real(const Vector3& vpos)> Sdf3::Box(const Vector3& vsize, const Vector3& vcenter) {
    fun = [vsize, vcenter](const Vector3& vpos) -> real {
        Vector3 vdist = (vpos - vcenter) - vsize * 0.5f;
        return (vdist.getMax(0)).len() + std::min(vdist.getMaxNum(), 0.f);
    };
    return fun;
}

std::function<real(const Vector3& vpos)> Sdf3::RoundBox(const Vector3& vsize, real fradius) {
    fun = [vsize, fradius](const Vector3 vpos) -> real {
        Vector3 vdist = vpos.getAbs() - vsize * 0.5f + fradius;
        return (vdist.getMax(0)).len() + std::fmin(vdist.getMaxNum(), 0) - fradius;
    };
    return fun;
}

std::function<real(const Vector3& vpos)> Sdf3::WireframeBox(const Vector3& vsize, real fthickness) {
    std::function<real(const Vector3 v)> funG = [&](Vector3 v) -> real {
        return Vector3(std::fmax(v.x, 0), std::fmax(v.y, 0), std::fmax(v.z, 0)).len() + std::fmin(std::fmax(v.x, std::fmax(v.y, v.z)), 0);
    };

    fun = [funG, vsize, fthickness](const Vector3 vpos) -> real {
        Vector3 vp = vpos.getAbs() - vsize / 2 - fthickness / 2;
        Vector3 vq = (vp + fthickness / 2).getAbs() - fthickness / 2;
        return std::fmin(std::fmin(funG(Vector3(vp.x, vq.y, vq.z)), funG(Vector3(vq.x, vp.y, vq.z))), funG(Vector3(vq.x, vq.y, vp.z)));
    };
    return fun;
}

std::function<real(const Vector3& vpos)> Sdf3::Torus(real fradius1, real fradius2) {
    fun = [fradius1, fradius2](const Vector3 vpos) -> real {
        real a = Vector2(vpos.x, vpos.y).len() - fradius1;
        return Vector2(a, vpos.z).len() - fradius2;
    };
    return fun;
}

std::function<real(const Vector3& vpos)> Sdf3::Capsule(const Vector3& vposa, const Vector3& vposb, real fradius) {
    fun = [vposa, vposb, fradius](const Vector3& vpos) -> real {
        Vector3 vpa = vpos - vposa;
        Vector3 vba = vposb - vposa;
        real h = Math::Clamp(vpa.dot(vba) / vba.dot(vba), 0.f, 1.f);
        return (vpa - vba * h).len() - fradius;
    };
    return fun;
}

std::function<real(const Vector3& vpos)> Sdf3::Cylinder(real fradius) {
    fun = [fradius](const Vector3& vpos) -> real {
        return Vector2(vpos.x, vpos.y).len() - fradius;
    };
    return fun;
}

std::function<real(const Vector3& vpos)> Sdf3::CappedCylinder(const Vector3& vposa, const Vector3& vposb, real fradius) {
    fun = [vposa, vposb, fradius](const Vector3& vpos) -> real {
        Vector3 vba = vposb - vposa;
        real baba = vba.dot(vba);
        real paba = vposa.dot(vba);
        real x = (vposa * baba - vba * paba).len() - fradius * baba;
        real y = std::fabs(paba - baba * 0.5f) - baba * 0.5f;
        real d = (std::fmax(x, y) < 0 ? -std::fmin(x * x, y * y * baba)
            : (x > 0 ? x * x : 0) + (y > 0 ? y * y * baba : 0));
        return std::signbit(d) * std::sqrt(std::fabs(d)) / baba;
    };
    return fun;
}

std::function<real(const Vector3& vpos)> Sdf3::RoundedCylinder(real fradiusa, real fradiusb, real height) {
    fun = [fradiusa, fradiusb, height](const Vector3& vpos) -> real {
        Vector2 vd = Vector2(Vector2(vpos.x, vpos.y).len() - fradiusa + fradiusb, std::fabs(vpos.z) - height / 2.f + fradiusb);
        return std::fmin(std::fmax(vd.x, vd.y), 0) + (Vector2(std::fmax(vd.x, 0), std::fmax(vd.y, 0)).len() - fradiusb);
    };
    return fun;
}

std::function<real(const Vector3& vpos)> Sdf3::CappedCone(const Vector3& vposa, const Vector3& vposb, real fradiusa, real fradiusb) {
    fun = [vposa, vposb, fradiusa, fradiusb](const Vector3& vpos) -> real {
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
    return fun;
}

std::function<real(const Vector3& vpos)> Sdf3::RoundedCone(real fradius1, real fradius2, real height) {
    fun = [fradius1, fradius2, height](const Vector3& vpos) -> real {
        Vector2 vq = Vector2(Vector2(vpos.x, vpos.y).len(), vpos.z);
        real fb = (fradius1 - fradius2) / height;
        real fa = std::sqrt(1 - fb * fb);
        real fk = vq.dot(Vector2(-fb, fa));
        real c1 = vq.len() - fradius1;
        real c2 = (vq - Vector2(0, height)).len() - fradius2;
        real c3 = vq.dot(Vector2(fa, fb)) - fradius1;
        return (fk < 0 ? c1 : (fk > fa * height ? c2 : c3));
    };
    return fun;
}

std::function<real(const Vector3& vpos)> Sdf3::Ellipsoid(real fsize) {
    fun = [fsize](const Vector3& vpos) -> real {
        real k0 = (vpos / fsize).len();
        real k1 = (vpos / (fsize * fsize)).len();
        return k0 * (k0 - 1) / k1;
    };
    return fun;
}

std::function<real(const Vector3& vpos)> Sdf3::Pyramid(real height) {
    fun = [height](const Vector3& vpos) -> real {
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
    return fun;
}

std::function<real(const Vector3& vpos)> Sdf3::Tetrahedron(real fradius) {
    fun = [fradius](const Vector3& vpos) -> real {
        return (std::max(std::fabs(vpos.x + vpos.y) - vpos.z, std::fabs(vpos.x - vpos.y) + vpos.z) - 1) / std::sqrt(3);
    };
    return fun;
}

std::function<real(const Vector3& vpos)> Sdf3::Octahedron(real fradius) {
    fun = [fradius](const Vector3& vpos) -> real {
        return (vpos.getAbs().getSumNum() - fradius) * std::tan(Math::PI / 60);
    };
    return fun;
}

std::function<real(const Vector3& vpos)> Sdf3::Dodecahedron(real fradius) {
    Vector3 v = Vector3((std::sqrt(5) + 1) / 2.f, 1.f, 0.f);
    fun = [v, fradius](const Vector3& vpos) -> real {
        Vector3 vp = (vpos / fradius).getAbs();
        real a = vp.dot(v);
        real b = vp.dot(Vector3(v.z, v.x, v.y));
        real c = vp.dot(Vector3(v.y, v.z, v.x));
        real q = (std::fmax(std::fmax(a, b), c) - v.x) * fradius;
        return q;
    };
    return fun;
}

std::function<real(const Vector3& vpos)> Sdf3::Icosahedron(real fradius) {
    Vector3 v = Vector3((std::sqrt(5) + 3) / 2.f, 1.f, 0.f);
    real w = std::sqrt(3) / 3;
    fun = [v, w, fradius](const Vector3& vpos) -> real {
        Vector3 vp = (vpos / fradius).getAbs();
        real a = vp.dot(v);
        real b = vp.dot(Vector3(v.z, v.x, v.y));
        real c = vp.dot(Vector3(v.y, v.z, v.x));
        real d = vp.dot(Vector3(w));
        return std::fmax(std::fmax(std::fmax(a, b), c) - v.x, d) * fradius;
    };
    return fun;
}

std::function<real(const Vector3& vpos)> Sdf3::TriPrism(Vector2 height) {
    fun = [height](const Vector3& vpos) -> real {
        Vector3 vq = vpos.getAbs();
        return std::fmax(vq.z - height.y, std::fmax(vq.x * 0.866025 + vpos.y * 0.5, -vpos.y) - height.x * 0.5f);
    };
    return fun;
}

// --------------------------------------------------------------
// operation
// --------------------------------------------------------------

std::function<real(const Vector3& vpos)> Sdf3::translate(std::function<real(const Vector3& vpos)> funSdf, const Vector3& voffset) {
    fun = [funSdf, voffset](const Vector3& vpos) -> real {
        return funSdf(vpos - voffset);
    };
    return fun;
}

std::function<real(const Vector3& vpos)> Sdf3::scale(std::function<real(const Vector3& vpos)> funSdf, const Vector3& vscale) {
    real m = std::min(vscale.x, std::min(vscale.y, vscale.z));
    fun = [&](const Vector3 vpos) -> real {
        return funSdf(vpos / vscale) * m;
    };
    return fun;
}

std::function<real(const Vector3& vpos)> Sdf3::rotate(std::function<real(const Vector3& vpos)> funSdf, real fangle, const Vector3& vdir) {
    Vector3 vdirn = vdir.getNormalize();
    real s = std::sin(fangle);
    real c = std::cos(fangle);
    real m = 1 - c;/*
    Eigen::Matrix3f matrix;
    matrix << m * vdir.x * vdir.x + c, m* vdir.x* vdir.y + vdir.z * s, m* vdir.z* vdir.x - vdir.y * s,
        m* vdir.x* vdir.y - vdir.z * s, m* vdir.y* vdir.y + c, m* vdir.y* vdir.z + vdir.x * s,
        m* vdir.z* vdir.x + vdir.y * s, m* vdir.y* vdir.z - vdir.x * s, m* vdir.z* vdir.z + c;
    matrix = matrix.transpose();
    */
    fun = [&](const Vector3 vpos) -> real {
        //Eigen::Vector3f vecpos(vpos.x, vpos.y, vpos.z);
        //return funSdf(Vector3(vecpos.dot(matrix.col(0)), vecpos.dot(matrix.col(1)), vecpos.dot(matrix.col(2))));
        return 0;
    };
    return fun;
}

std::function<real(const Vector3& vpos)> Sdf3::rotate_to(std::function<real(const Vector3& vpos)> funSdf, const Vector3& va, const Vector3& vb) {
    real angle = va.dot(vb);
    Vector3 vdir = vb.cross(va);
    return rotate(funSdf, angle, vdir);
}

std::function<real(const Vector3& vpos)> Sdf3::orient(std::function<real(const Vector3& vpos)> funSdf, const Vector3& vaxis) {
    return rotate_to(funSdf, Vector3(0, 1, 0), vaxis);
}

std::function<real(const Vector3& vpos)> Sdf3::circular_array(std::function<real(const Vector3& vpos)> funSdf, int count) {
    int da = 2 * Math::PI / count;
    fun = [&](const Vector3 vpos) -> real {
        real d = Vector2(vpos.x, vpos.y).len();
        real a = std::atan(vpos.y / vpos.x);
        real d1 = funSdf(Vector3(a - da, std::sin(a - da), vpos.z));
        real d2 = funSdf(Vector3(std::cos(a) * d, std::sin(a) * d, vpos.z));
        return std::min(d1, d2);
    };
    return fun;
}

std::function<real(const Vector3& vpos)> Sdf3::elongate(std::function<real(const Vector3& vpos)> funSdf, const Vector3& vsize) {
    fun = [&](const Vector3 vpos) -> real {
        Vector3 vq = vpos.getAbs() - vsize;
        real w = std::min(std::max(vq.x, std::max(vq.y, vq.z)), 0.0f);
        return funSdf(vq.getMax(0)) + w;
    };
    return fun;
}

std::function<real(const Vector3& vpos)> Sdf3::twist(std::function<real(const Vector3& vpos)> funSdf, int k) {
    fun = [&](const Vector3 vpos) -> real {
        real c = std::cos(k * vpos.z);
        real s = std::sin(k * vpos.z);
        Vector3 vres = Vector3(c * vpos.x - s * vpos.y, s * vpos.x + c * vpos.y, vpos.z);
        return funSdf(vres);
    };
    return fun;
}

std::function<real(const Vector3& vpos)> Sdf3::bend(std::function<real(const Vector3& vpos)> funSdf, int k) {
    fun = [&](const Vector3& vpos) -> real {
        real c = std::cos(k * vpos.x);
        real s = std::sin(k * vpos.x);
        Vector3 vres = Vector3(c * vpos.x - s * vpos.y, s * vpos.x + c * vpos.y, vpos.z);
        return funSdf(vres);
    };
    return fun;
}

std::function<real(const Vector3& vpos)> Sdf3::bend_linear(std::function<real(const Vector3& vpos)> funSdf, const Vector3& vpa, const Vector3& vpb, const Vector3& v) {
    Vector3 vab = (vpb - vpa);
    fun = [&](const Vector3 vpos) -> real {
        real t = Math::Clamp((vpos - vpa).dot(vab) / vab.dot(vab), 0.f, 1.f);
        return funSdf(vpos + t * v);
    };
    return fun;
}

std::function<real(const Vector3& vpos)> Sdf3::bend_radial(std::function<real(const Vector3& vpos)> funSdf) {
    fun = [&](const Vector3 vpos) -> real {
        return 0;
    };
    return fun;
}
std::function<real(const Vector3& vpos)> Sdf3::transition_linear(std::function<real(const Vector3& vpos)> funSdf) {
    fun = [&](const Vector3 vpos) -> real {
        return 0;
    };
    return fun;
}
std::function<real(const Vector3& vpos)> Sdf3::transition_radial(std::function<real(const Vector3& vpos)> funSdf) {
    fun = [&](const Vector3 vpos) -> real {
        return 0;
    };
    return fun;
}
std::function<real(const Vector3& vpos)> Sdf3::wrap_around(std::function<real(const Vector3& vpos)> funSdf) {
    fun = [&](const Vector3 vpos) -> real {
        return 0;
    };
    return fun;
}

std::function<real(const Vector3& vpos)> Sdf3::slice(std::function<real(const Vector3& vpos)> funSdf) {
    Slab(Vector3(0, 0, -1), Vector3(0, 0, 1));
    fun = [&](const Vector3 vpos) -> real {
        return 0;
    };
    return fun;
}

// --------------------------------------------------------------
// dn
// --------------------------------------------------------------

std::function<real(const Vector3& vpos)> Sdf3::Union(vector<std::function<real(const Vector3& vpos)> > funs) {
    std::function<real(const Vector3& vpos)> fun = [&](const Vector3 vpos) -> real {
        return 0;
    };
    return fun;
}
std::function<real(const Vector3& vpos)> Sdf3::difference(vector<std::function<real(const Vector3& vpos)> > funs) {
    fun = [&](const Vector3& vpos) -> real {
        return 0;
    };
    return fun;
}
std::function<real(const Vector3& vpos)> Sdf3::intersection(vector<std::function<real(const Vector3& vpos)> > funs, float k) {
    std::function<real(const Vector3 vpos)> fun = [&](const Vector3 vpos) -> real {
        real d1 = 0;
        for (auto f : funs) {
            real d2 = f(vpos);
            float h = Math::Clamp(0.5f - 0.5f * (d2 - d1) / k, 0.f, 1.f);
            real m = d2 + (d1 - d2) * h;
            d1 = m + k * h * (1 - h);
        }
        return d1;
    };
    return fun;
}

std::function<real(const Vector3& vpos)> Sdf3::blend(vector<std::function<real(const Vector3& vpos)> > funs) {
    std::function<real(const Vector3 vpos)> fun = [&](const Vector3 vpos) -> real {
        return 0;
    };
    return fun;
}
std::function<real(const Vector3& vpos)> Sdf3::negate(vector<std::function<real(const Vector3& vpos)> > funs) {
    std::function<real(const Vector3 vpos)> fun = [&](const Vector3 vpos) -> real {
        return 0;
    };
    return fun;
}
std::function<real(const Vector3& vpos)> Sdf3::dilate(vector<std::function<real(const Vector3& vpos)> > funs) {
    std::function<real(const Vector3 vpos)> fun = [&](const Vector3 vpos) -> real {
        return 0;
    };
    return fun;
}
std::function<real(const Vector3& vpos)> Sdf3::erode(std::function<real(const Vector3& vpos)> funSdf, real r) {
    std::function<real(const Vector3 vpos)> fun = [&](const Vector3 vpos) -> real {
        return 0;
    };
    return fun;
}
std::function<real(const Vector3& vpos)> Sdf3::shell(std::function<real(const Vector3& vpos)> funSdf, real r) {
    std::function<real(const Vector3 vpos)> fun = [&](const Vector3 vpos) -> real {
        return 0;
    };
    return fun;
}
std::function<real(const Vector3& vpos)> Sdf3::repeat(std::function<real(const Vector3& vpos)> funSdf, real spacing) {
    std::function<real(const Vector3 vpos)> fun = [&](const Vector3 vpos) -> real {
        return 0;
    };
    return fun;
}
