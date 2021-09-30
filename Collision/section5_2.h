#ifndef SECTION5_2_H
#define SECTION5_2_H

#include "section5.h"



//----------------------------------------
//-------------------- section 5.2 --------------------
//----------------------------------------



// Determine whether plane p intersects sphere s
int TestSpherePlane(Sphere s, Plane p)
{
    // For a normalized plane (|p.n| = 1), evaluating the plane equation
    // for a point gives the signed distance of the point to the plane
    float dist = Dot(s.c, p.n) - p.d;
    // If sphere center within +/-radius from plane, plane intersects sphere
    return Abs(dist) <= s.r;
}

// Determine whether sphere s fully behind (inside negative halfspace of) plane p
int InsideSpherePlane(Sphere s, Plane p)
{
    float dist = Dot(s.c, p.n) - p.d;
    return dist < -s.r;
}

// Determine whether sphere s intersects negative halfspace of plane p
int TestSphereHalfspace(Sphere s, Plane p)
{
    float dist = Dot(s.c, p.n) - p.d;
    return dist <= s.r;
}

// Test if OBB b intersects plane p
int TestOBBPlane(OBB b, Plane p)
{
    // Compute the projection interval radius of b onto L(t) = b.c + t * p.n
    float r = b.e[0]*Abs(Dot(p.n, b.u[0])) +
              b.e[1]*Abs(Dot(p.n, b.u[1])) +
              b.e[2]*Abs(Dot(p.n, b.u[2]));
    // Compute distance of box center from plane
    float s = Dot(p.n, b.c) - p.d;
    // Intersection occurs when distance s falls within [-r,+r] interval
    return Abs(s) <= r;
}

// Test if AABB b intersects plane p
int TestAABBPlane(AABB b, Plane p)
{
    // These two lines not necessary with a (center, extents) AABB representation
    Vector3 c = (b.max + b.min) * 0.5f; // Compute AABB center
    Vector3 e = b.max - c; // Compute positive extents

    // Compute the projection interval radius of b onto L(t) = b.c + t * p.n
    float r = e[0]*Abs(p.n[0]) + e[1]*Abs(p.n[1]) + e[2]*Abs(p.n[2]);
    // Compute distance of box center from plane
    float s = Dot(p.n, c) - p.d;
    // Intersection occurs when distance s falls within [-r,+r] interval
    return Abs(s) <= r;
}

// Returns true if sphere s intersects AABB b, false otherwise
int TestSphereAABB(Sphere s, AABB b)
{
    // Compute squared distance between sphere center and AABB
    float sqDist = SqDistPointAABB(s.c, b);

    // Sphere and AABB intersect if the (squared) distance
    // between them is less than the (squared) sphere radius
    return sqDist <= s.r * s.r;
}

// Returns true if sphere s intersects AABB b, false otherwise.
// The point p on the AABB closest to the sphere center is also returned
int TestSphereAABB(Sphere s, AABB b, Vector3 &p)
{
    // Find point p on AABB closest to sphere center
    ClosestPtPointAABB(s.c, b, p);

    // Sphere and AABB intersect if the (squared) distance from sphere
    // center to point p is less than the (squared) sphere radius
    Vector3 v = p - s.c;
    return Dot(v, v) <= s.r * s.r;
}

// Returns true if sphere s intersects OBB b, false otherwise.
// The point p on the OBB closest to the sphere center is also returned
int TestSphereOBB(Sphere s, OBB b, Vector3 &p)
{
    // Find point p on OBB closest to sphere center
    ClosestPtPointOBB(s.c, b, p);

    // Sphere and OBB intersect if the (squared) distance from sphere
    // center to point p is less than the (squared) sphere radius
    Vector3 v = p - s.c;
    return Dot(v, v) <= s.r * s.r;
}

// Returns true if sphere s intersects triangle ABC, false otherwise.
// The point p on abc closest to the sphere center is also returned
int TestSphereTriangle(Sphere s, Vector3 a, Vector3 b, Vector3 c, Vector3 &p)
{
    // Find point P on triangle ABC closest to sphere center
    p = ClosestPtPointTriangle(s.c, a, b, c);

    // Sphere and triangle intersect if the (squared) distance from sphere
    // center to point p is less than the (squared) sphere radius
    Vector3 v = p - s.c;
    return Dot(v, v) <= s.r * s.r;
}

//function not define
/*
struct Polygon {
    std::vector<Vector3> v;
    int numVerts;
};

// Test whether sphere s intersects polygon p
int TestSpherePolygon(Sphere s, Polygon p)
{
    // Compute normal for the plane of the polygon
    Vector3 n = Normalize(Cross(p.v[1] - p.v[0], p.v[2] - p.v[0]));
    // Compute the plane equation for p
    Plane m; m.n = n; m.d = -Dot(n, p.v[0]);
    // No intersection if sphere not intersecting plane of polygon
    if (!TestSpherePlane(s, m)) return 0;
    // Test to see if any one of the polygon edges pierces the sphere
    for (int k = p.numVerts, i = 0, j = k - 1; i < k; j = i, i++) {
        float t;
        Vector3 q;
        // Test if edge (p.v[j], p.v[i]) intersects s
        if (IntersectRaySphere(p.v[j], p.v[i] - p.v[j], s, t, q) && t <= 1.0f)
            return 1;
    }
    // Test if the orthogonal projection q of the sphere center onto m is inside p
    Vector3 q = ClosestPtPointPlane(s.c, m);
    return PointInPolygon(q, p);
}
*/

//section 13
int TestTriangleAABB(Vector3  v0, Vector3  v1, Vector3  v2, AABB b)
{
    float p0, p1, p2, r;

    // Compute box center and extents (if not already given in that format)
    Vector3 c = (b.min + b.max) * 0.5f;
    float e0 = (b.max.x - b.min.x) * 0.5f;
    float e1 = (b.max.y - b.min.y) * 0.5f;
    float e2 = (b.max.z - b.min.z) * 0.5f;

    // Translate triangle as conceptually moving AABB to origin
    v0 = v0 - c;
    v1 = v1 - c;
    v2 = v2 - c;

    // Compute edge vectors for triangle
    Vector3 f0 = v1 - v0,  f1 = v2 - v1, f2 = v0 - v2;

    // Test axes a00..a22 (category 3)
    // Test axis a00
    p0 = v0.z*v1.y - v0.y*v1.z;
    p2 = v2.z*(v1.y - v0.y) - v2.y*(v1.z - v0.z);
    r = e1 * Abs(f0.z) + e2 * Abs(f0.y);
    if (Max(-Max(p0, p2), Min(p0, p2)) > r) return 0; // Axis is a separating axis

    // Repeat similar tests for remaining axes a01..a22
    //...

    // Test the three axes corresponding to the face normals of AABB b (category 1).
    // Exit if...
    // ... [-e0, e0] and [min(v0.x,v1.x,v2.x), max(v0.x,v1.x,v2.x)] do not overlap
    if (Max(v0.x, v1.x, v2.x) < -e0 || Min(v0.x, v1.x, v2.x) > e0) return 0;
    // ... [-e1, e1] and [min(v0.y,v1.y,v2.y), max(v0.y,v1.y,v2.y)] do not overlap
    if (Max(v0.y, v1.y, v2.y) < -e1 || Min(v0.y, v1.y, v2.y) > e1) return 0;
    // ... [-e2, e2] and [min(v0.z,v1.z,v2.z), max(v0.z,v1.z,v2.z)] do not overlap
    if (Max(v0.z, v1.z, v2.z) < -e2 || Min(v0.z, v1.z, v2.z) > e2) return 0;

    // Test separating axis corresponding to triangle face normal (category 2)
    Plane p;
    p.n = Cross(f0, f1);
    p.d = Dot(p.n, v0);
    return TestAABBPlane(b, p);
}



//----------------------------------------
//-------------------- section 5.3 --------------------
//----------------------------------------



int IntersectSegmentPlane(Vector3 a, Vector3 b, Plane p, float &t, Vector3 &q)
{
    // Compute the t value for the directed line ab intersecting the plane
    Vector3 ab = b - a;
    t = (p.d - Dot(p.n, a)) / Dot(p.n, ab);

    // If t in [0..1] compute and return intersection point
    if (t >= 0.0f && t <= 1.0f) {
        q = a + t * ab;
        return 1;
    }
    // Else no intersection
    return 0;
}

// Intersect segment ab against plane of triangle def. If intersecting,
// return t value and position q of intersection
int IntersectSegmentPlane(Vector3 a, Vector3 b, Vector3 d, Vector3 e, Vector3 f,
                          float &t, Vector3 &q)
{
    Plane p;
    p.n = Cross(e - d, f - d);
    p.d = Dot(p.n, d);
    return IntersectSegmentPlane(a, b, p, t, q);
}

// Intersects ray r = p + td, |d| = 1, with sphere s and, if intersecting,
// returns t value of intersection and intersection point q
int IntersectRaySphere(Vector3 p, Vector3 d, Sphere s, float &t, Vector3 &q)
{
    Vector3 m = p - s.c;
    float b = Dot(m, d);
    float c = Dot(m, m) - s.r * s.r;
    // Exit if r's origin outside s (c > 0)and r pointing away from s (b > 0)
    if (c > 0.0f && b > 0.0f) return 0;
    float discr = b*b - c;
    // A negative discriminant corresponds to ray missing sphere
    if (discr < 0.0f) return 0;
    // Ray now found to intersect sphere, compute smallest t value of intersection
    t = -b - std::sqrt(discr);
    // If t is negative, ray started inside sphere so clamp t to zero
    if (t < 0.0f) t = 0.0f;
    q = p + t * d;
    return 1;
}

// Test if ray r = p + td intersects sphere s
int TestRaySphere(Vector3 p, Vector3 d, Sphere s)
{
    Vector3 m = p - s.c;
    float c = Dot(m, m) - s.r * s.r;
    // If there is definitely at least one real root, there must be an intersection
    if (c <= 0.0f) return 1;
    float b = Dot(m, d);
    // Early exit if ray origin outside sphere and ray pointing away from sphere
    if (b > 0.0f) return 0;
    float disc = b*b - c;
    // A negative discriminant corresponds to ray missing sphere
    if (disc < 0.0f) return 0;
    // Now ray must hit sphere
    return 1;
}

// Intersect ray R(t) = p + t*d against AABB a. When intersecting,
// return intersection distance tmin and point q of intersection
int IntersectRayAABB(Vector3 p, Vector3 d, AABB a, float &tmin, Vector3 &q)
{
    tmin = 0.0f;          // set to -FLT_MAX to get first hit on line
    float tmax = POS_INFINITY; // set to max distance ray can travel (for segment)

    // For all three slabs
    for (int i = 0; i < 3; i++) {
        if (Abs(d[i]) < EPSILON) {
            // Ray is parallel to slab. No hit if origin not within slab
            if (p[i] < a.min[i] || p[i] > a.max[i]) return 0;
        } else {
            // Compute intersection t value of ray with near and far plane of slab
            float ood = 1.0f / d[i];
            float t1 = (a.min[i] - p[i]) * ood;
            float t2 = (a.max[i] - p[i]) * ood;
            // Make t1 be intersection with near plane, t2 with far plane
            if (t1 > t2) Swap(t1, t2);
            // Compute the intersection of slab intersections intervals
            tmin = Max(tmin, t1);
            tmax = Min(tmax, t2);
            // Exit with no collision as soon as slab intersection becomes empty
            if (tmin > tmax) return 0;
        }
    }
    // Ray intersects all 3 slabs. Return point (q) and intersection t value (tmin)
    q = p + d * tmin;
    return 1;
}

// Test if segment specified by points p0 and p1 intersects AABB b
int TestSegmentAABB(Vector3 p0, Vector3 p1, AABB b)
{
    Vector3 c = (b.min + b.max) * 0.5f; // Box center-point
    Vector3 e = b.max - c;             // Box halflength extents
    Vector3 m = (p0 + p1) * 0.5f;       // Segment midpoint
    Vector3 d = p1 - m;                // Segment halflength vector
    m = m - c;                        // Translate box and segment to origin

    // Try world coordinate axes as separating axes
    float adx = Abs(d.x);
    if (Abs(m.x) > e.x + adx) return 0;
    float ady = Abs(d.y);
    if (Abs(m.y) > e.y + ady) return 0;
    float adz = Abs(d.z);
    if (Abs(m.z) > e.z + adz) return 0;

    // Add in an epsilon term to counteract arithmetic errors when segment is
    // (near) parallel to a coordinate axis (see text for detail)
    adx += EPSILON; ady += EPSILON; adz += EPSILON;

    // Try cross products of segment direction vector with coordinate axes
    if (Abs(m.y * d.z - m.z * d.y) > e.y * adz + e.z * ady) return 0;
    if (Abs(m.z * d.x - m.x * d.z) > e.x * adz + e.z * adx) return 0;
    if (Abs(m.x * d.y - m.y * d.x) > e.x * ady + e.y * adx) return 0;

    // No separating axis found; segment must be overlapping AABB
    return 1;
}

float ScalarTriple(const Vector3 &a, const Vector3 &b, const Vector3 &c) {
    return a.cross(b).dot(c);
}

// Given line pq and ccw triangle abc, return whether line pierces triangle. If
// so, also return the barycentric coordinates (u,v,w) of the intersection point
int IntersectLineTriangle(Vector3 p, Vector3 q, Vector3 a, Vector3 b, Vector3 c,
                          float &u, float &v, float &w)
{
    Vector3 pq = q - p;
    Vector3 pa = a - p;
    Vector3 pb = b - p;
    Vector3 pc = c - p;
    // Test if pq is inside the edges bc, ca and ab. Done by testing
    // that the signed tetrahedral volumes, computed using scalar triple
    // products, are all positive
    u = ScalarTriple(pq, pc, pb);
    if (u < 0.0f) return 0;
    v = ScalarTriple(pq, pa, pc);
    if (v < 0.0f) return 0;
    w = ScalarTriple(pq, pb, pa);
    if (w < 0.0f) return 0;

    // Compute the barycentric coordinates (u, v, w) determining the
    // intersection point r, r = u*a + v*b + w*c
    float denom = 1.0f / (u + v + w);
    u *= denom;
    v *= denom;
    w *= denom; // w = 1.0f - u - v;
    return 1;
}

// Given line pq and ccw quadrilateral abcd, return whether the line
// pierces the triangle. If so, also return the point r of intersection
int IntersectLineQuad(Vector3 p, Vector3 q, Vector3 a, Vector3 b, Vector3 c, Vector3 d, Vector3 &r)
{
    Vector3 pq = q - p;
    Vector3 pa = a - p;
    Vector3 pb = b - p;
    Vector3 pc = c - p;
    // Determine which triangle to test against by testing against diagonal first
    Vector3 m = Cross(pc, pq);
    float v = Dot(pa, m); // ScalarTriple(pq, pa, pc);
    if (v >= 0.0f) {
        // Test intersection against triangle abc
        float u = -Dot(pb, m); // ScalarTriple(pq, pc, pb);
        if (u < 0.0f) return 0;
        float w = ScalarTriple(pq, pb, pa);
        if (w < 0.0f) return 0;
        // Compute r, r = u*a + v*b + w*c, from barycentric coordinates (u, v, w)
        float denom = 1.0f / (u + v + w);
        u *= denom;
        v *= denom;
        w *= denom; // w = 1.0f - u - v;
        r = u*a + v*b + w*c;
    } else {
        // Test intersection against triangle dac
        Vector3 pd = d - p;
        float u = Dot(pd, m); // ScalarTriple(pq, pd, pc);
        if (u < 0.0f) return 0;
        float w = ScalarTriple(pq, pa, pd);
        if (w < 0.0f) return 0;
        v = -v;
        // Compute r, r = u*a + v*d + w*c, from barycentric coordinates (u, v, w)
        float denom = 1.0f / (u + v + w);
        u *= denom;
        v *= denom;
        w *= denom; // w = 1.0f - u - v;
        r = u*a + v*d + w*c;
    }
    return 1;
}

// Given segment pq and triangle abc, returns whether segment intersects
// triangle and if so, also returns the barycentric coordinates (u,v,w)
// of the intersection point
int IntersectSegmentTriangle(Vector3 p, Vector3 q, Vector3 a, Vector3 b, Vector3 c,
                             float &u, float &v, float &w, float &t)
{
    Vector3 ab = b - a;
    Vector3 ac = c - a;
    Vector3 qp = p - q;

    // Compute triangle normal. Can be precalculated or cached if
    // intersecting multiple segments against the same triangle
    Vector3 n = Cross(ab, ac);

    // Compute denominator d. If d <= 0, segment is parallel to or points
    // away from triangle, so exit early
    float d = Dot(qp, n);
    if (d <= 0.0f) return 0;

    // Compute intersection t value of pq with plane of triangle. A ray
    // intersects iff 0 <= t. Segment intersects iff 0 <= t <= 1. Delay
    // dividing by d until intersection has been found to pierce triangle
    Vector3 ap = p - a;
    t = Dot(ap, n);
    if (t < 0.0f) return 0;
    if (t > d) return 0; // For segment; exclude this code line for a ray test

    // Compute barycentric coordinate components and test if within bounds
    Vector3 e = Cross(qp, ap);
    v = Dot(ac, e);
    if (v < 0.0f || v > d) return 0;
    w = -Dot(ab, e);
    if (w < 0.0f || v + w > d) return 0;

    // Segment/ray intersects triangle. Perform delayed division and
    // compute the last barycentric coordinate component
    float ood = 1.0f / d;
    t *= ood;
    v *= ood;
    w *= ood;
    u = 1.0f - v - w;
    return 1;
}

struct Triangle {
    Plane p;           // Plane equation for triangle plane
    Plane edgePlaneBC; // When evaluated gives barycentric weight u (for vertex A)
    Plane edgePlaneCA; // When evaluated gives barycentric weight v (for vertex B)
};

// Given segment pq and precomputed triangle tri, returns whether segment intersects
// triangle. If so, also returns the barycentric coordinates (u,v,w) of the
// intersection point s, and the parameterized intersection t value
int IntersectSegmentTriangle(Vector3 p, Vector3 q, Triangle tri,
                             float &u, float &v, float &w, float &t, Vector3 &s)
{
    // Compute distance of p to triangle plane. Exit if p lies behind plane
    float distp = Dot(p, tri.p.n) - tri.p.d;
    if (distp < 0.0f) return 0;

    // Compute distance of q to triangle plane. Exit if q lies in front of plane
    float distq = Dot(q, tri.p.n) - tri.p.d;
    if (distq >= 0.0f) return 0;

    // Compute t value and point s of intersection with triangle plane
    float denom = distp - distq;
    t = distp / denom;
    s = p + t * (q - p);

    // Compute the barycentric coordinate u; exit if outside 0..1 range
    u = Dot(s, tri.edgePlaneBC.n) - tri.edgePlaneBC.d;
    if (u < 0.0f || u > 1.0f) return 0;
    // Compute the barycentric coordinate v; exit if negative
    v = Dot(s, tri.edgePlaneCA.n) - tri.edgePlaneCA.d;
    if (v < 0.0f) return 0;
    // Compute the barycentric coordinate w; exit if negative
    w = 1.0f - u - v;
    if (w < 0.0f) return 0;

    // Segment intersects tri at distance t in position s (s = u*A + v*B + w*C)
    return 1;
}

// Intersect segment S(t)=sa+t(sb-sa), 0<=t<=1 against cylinder specified by p, q and r
int IntersectSegmentCylinder(Vector3 sa, Vector3 sb, Vector3 p, Vector3 q, float r, float &t)
{
    Vector3 d = q - p, m = sa - p, n = sb - sa;
    float md = Dot(m, d);
    float nd = Dot(n, d);
    float dd = Dot(d, d);
    // Test if segment fully outside either endcap of cylinder
    if (md < 0.0f && md + nd < 0.0f) return 0; // Segment outside 憄- side of cylinder
    if (md > dd && md + nd > dd) return 0;     // Segment outside 憅- side of cylinder
    float nn = Dot(n, n);
    float mn = Dot(m, n);
    float a = dd * nn - nd * nd;
    float k = Dot(m, m) - r * r;
    float c = dd * k - md * md;
    if (Abs(a) < EPSILON) {
        // Segment runs parallel to cylinder axis
        if (c > 0.0f) return 0; // 慳- and thus the segment lie outside cylinder
        // Now known that segment intersects cylinder; figure out how it intersects
        if (md < 0.0f) t = -mn / nn; // Intersect segment against 憄- endcap
        else if (md > dd) t = (nd - mn) / nn; // Intersect segment against 憅- endcap
        else t = 0.0f; // 慳- lies inside cylinder
        return 1;
    }
    float b = dd * mn - nd * md;
    float discr = b * b - a * c;
    if (discr < 0.0f) return 0; // No real roots; no intersection
    t = (-b - std::sqrt(discr)) / a;
    if (t < 0.0f || t > 1.0f) return 0; // Intersection lies outside segment
    if (md + t * nd < 0.0f) {
        // Intersection outside cylinder on 憄- side
        if (nd <= 0.0f) return 0; // Segment pointing away from endcap
        t = -md / nd;
        // Keep intersection if Dot(S(t) - p, S(t) - p) <= r^2
        return k + 2 * t * (mn + t * nn) <= 0.0f;
    } else if (md + t * nd > dd) {
        // Intersection outside cylinder on 憅- side
        if (nd >= 0.0f) return 0; // Segment pointing away from endcap
        t = (dd - md) / nd;
        // Keep intersection if Dot(S(t) - q, S(t) - q) <= r^2
        return k + dd - 2 * md + t * (2 * (mn - nd) + t * nn) <= 0.0f;
    }
    // Segment intersects cylinder between the end-caps; t is correct
    return 1;
}

// Intersect segment S(t)=A+t(B-A), 0<=t<=1 against convex polyhedron specified
// by the n halfspaces defined by the planes p[]. On exit tfirst and tlast
// define the intersection, if any
int IntersectSegmentPolyhedron(Vector3 a, Vector3 b, Plane p[], int n,
                               float &tfirst, float &tlast)
{
    // Compute direction vector for the segment
    Vector3 d = b - a;
    // Set initial interval to being the whole segment. For a ray, tlast should be
    // set to +FLT_MAX. For a line, additionally tfirst should be set to 朏LT_MAX
    tfirst = 0.0f;
    tlast = 1.0f;
    // Intersect segment against each plane
    for (int i = 0; i < n; i++) {
        float denom = Dot(p[i].n, d);
        float dist = p[i].d - Dot(p[i].n, a);
        // Test if segment runs parallel to the plane
        if (denom == 0.0f) {
            // If so, return 搉o intersection- if segment lies outside plane
            if (dist > 0.0f) return 0;
        } else {
            // Compute parameterized t value for intersection with current plane
            float t = dist / denom;
            if (denom < 0.0f) {
                // When entering halfspace, update tfirst if t is larger
                if (t > tfirst) tfirst = t;
            } else {
                // When exiting halfspace, update tlast if t is smaller
                if (t < tlast) tlast = t;
            }
            // Exit with 搉o intersection- if intersection becomes empty
            if (tfirst > tlast) return 0;
        }
    }
    // A nonzero logical intersection, so the segment intersects the polyhedron
    return 1;
}

#endif // SECTION5_2_H
