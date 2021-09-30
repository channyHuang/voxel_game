#ifndef SECTION5_3_H
#define SECTION5_3_H

#include "section5.h"


//----------------------------------------
//-------------------- section 5.4 --------------------
//----------------------------------------



inline float get_normal_angle(const Vector3 &normal1, const Vector3 &normal2) {
    if (normal1.len() == 0 || normal2.len() == 0) return 0.f;
    float cos_theta = normal1.dot(normal2); // / (normal1.len() * normal2.len());
    return acos(cos_theta); //[0,pi]
}
//2d point a - b - c is anti cross wise
//cross > 0 - anti cross wise
static int TriangleIsCCW(const Vector2 &a, const Vector2 &b, const Vector2 &c) {
    float cross = ((a.x - c.x) * (b.y - c.y)) - ((a.y - c.y) * (b.x - c.x));
    return (cross == 0 ? 0 : (cross > 0 ? 1 : -1));
}
bool TriangleIsCCW(const Vector3 &a, const Vector3 &b, const Vector3 &c) {
// > 0 for counterclockwise
// = 0 for none (degenerate)
// < 0 for clockwise
// Triangle vertices
    Vector3 n1 = (b - a).cross(c - a);
    Vector3 n2 = (c - b).cross(a - b);
    return get_normal_angle(n1, n2) < PI / 2;  // Small epsilon value
}

// Test if point p lies inside ccw-specified convex n-gon given by vertices v[]
int PointInConvexPolygon(Vector3 p, int n, Vector3 v[])
{
    // Do binary search over polygon vertices to find the fan triangle
    // (v[0], v[low], v[high]) the point p lies within the near sides of
    int low = 0, high = n;
    do {
        int mid = (low + high) / 2;
        if (TriangleIsCCW(v[0], v[mid], p))
            low = mid;
        else
            high = mid;
    } while (low + 1 < high);

    // If point outside last (or first) edge, then it is not inside the n-gon
    if (low == 0 || high == n) return 0;

    // p is inside the polygon if it is left of
    // the directed edge from v[low] to v[high]
    return TriangleIsCCW(v[low], v[high], p);
}

// Test if point P lies inside the counterclockwise triangle ABC
int PointInTriangle(Vector3 p, Vector3 a, Vector3 b, Vector3 c)
{
    // Translate point and triangle so that point lies at origin
    a -= p; b -= p; c -= p;
    // Compute normal vectors for triangles pab and pbc
    Vector3 u = Cross(b, c);
    Vector3 v = Cross(c, a);
    // Make sure they are both pointing in the same direction
    if (Dot(u, v) < 0.0f) return 0;
    // Compute normal vector for triangle pca
    Vector3 w = Cross(a, b);
    // Make sure it points in the same direction as the first two
    if (Dot(u, w) < 0.0f) return 0;
    // Otherwise P must be in (or on) the triangle
    return 1;
}

namespace Section5 {
    // Test if point P lies inside the counterclockwise 3D triangle ABC
    int PointInTriangle(Vector3 p, Vector3 a, Vector3 b, Vector3 c)
    {
        // Translate point and triangle so that point lies at origin
        a -= p; b -= p; c -= p;

        float ab = Dot(a, b);
        float ac = Dot(a, c);
        float bc = Dot(b, c);
        float cc = Dot(c, c);
        // Make sure plane normals for pab and pbc point in the same direction
        if (bc * ac - cc * ab < 0.0f) return 0;
        // Make sure plane normals for pab and pca point in the same direction
        float bb = Dot(b, b);
        if (ab * bc - ac * bb < 0.0f) return 0;
        // Otherwise P must be in (or on) the triangle
        return 1;
    }
}

// Compute the 2D pseudo cross product Dot(Perp(u), v)
float Cross2D(Vector2 u, Vector2 v)
{
    return u.y * v.x - u.x * v.y;
}

// Test if 2D point P lies inside the counterclockwise 2D triangle ABC
int PointInTriangle(Vector2 p, Vector2 a, Vector2 b, Vector2 c)
{
    // If P to the right of AB then outside triangle
    if (Cross2D(p - a, b - a) < 0.0f) return 0;
    // If P to the right of BC then outside triangle
    if (Cross2D(p - b, c - b) < 0.0f) return 0;
    // If P to the right of CA then outside triangle
    if (Cross2D(p - c, a - c) < 0.0f) return 0;
    // Otherwise P must be in (or on) the triangle
    return 1;
}

bool SameSign(float a, float b) {
    return (a >= 0 && b >= 0) || (a <= 0 && b <= 0);
}

// Test if 2D point P lies inside 2D triangle ABC
int PointInTriangle2D(Vector2 p, Vector2 a, Vector2 b, Vector2 c)
{
    float pab = Cross2D(p - a, b - a);
    float pbc = Cross2D(p - b, c - b);
    // If P left of one of AB and BC and right of the other, not inside triangle
    if (!SameSign(pab, pbc)) return 0;
    float pca = Cross2D(p - c, a - c);
    // If P left of one of AB and CA and right of the other, not inside triangle
    if (!SameSign(pab, pca)) return 0;
    // P left or right of all edges, so must be in (or on) the triangle
    return 1;
}

// Test if point p inside polyhedron given as the intersection volume of n halfspaces
int TestPointPolyhedron(Vector3 p, Plane *h, int n)
{
    for (int i = 0; i < n; i++) {
        // Exit with 憂o containment- if p ever found outside a halfspace
        if (DistPointPlane(p, h[i]) > 0.0f) return 0;
    }
    // p inside all halfspaces, so p must be inside intersection volume
    return 1;
}

// Given planes p1 and p2, compute line L = p+t*d of their intersection.
// Return 0 if no such line exists
int IntersectPlanes(Plane p1, Plane p2, Vector3 &p, Vector3 &d)
{
    // Compute direction of intersection line
    d = Cross(p1.n, p2.n);

    // If d is zero, the planes are parallel (and separated)
    // or coincident, so they抮e not considered intersecting
    if (Dot(d, d) < EPSILON) return 0;

    float d11 = Dot(p1.n, p1.n);
    float d12 = Dot(p1.n, p2.n);
    float d22 = Dot(p2.n, p2.n);

    float denom = d11*d22 - d12*d12;
    float k1 = (p1.d*d22 - p2.d*d12) / denom;
    float k2 = (p2.d*d11 - p1.d*d12) / denom;
    p = k1*p1.n + k2*p2.n;
    return 1;
}

namespace Section5 {
    // Given planes p1 and p2, compute line L = p+t*d of their intersection.
    // Return 0 if no such line exists
    int IntersectPlanes(Plane p1, Plane p2, Vector3 &p, Vector3 &d)
    {
        // Compute direction of intersection line
        d = Cross(p1.n, p2.n);

        // If d is (near) zero, the planes are parallel (and separated)
        // or coincident, so they抮e not considered intersecting
        float denom = Dot(d, d);
        if (denom < EPSILON) return 0;

        // Compute point on intersection line
        p = Cross(p1.d*p2.n - p2.d*p1.n, d) / denom;
        return 1;
    }
}

// Compute the point p at which the three planes p1, p2 and p3 intersect (if at all)
int IntersectPlanes(Plane p1, Plane p2, Plane p3, Vector3 &p)
{
    Vector3 m1 = Vector3(p1.n.x, p2.n.x, p3.n.x);
    Vector3 m2 = Vector3(p1.n.y, p2.n.y, p3.n.y);
    Vector3 m3 = Vector3(p1.n.z, p2.n.z, p3.n.z);

    Vector3 u = Cross(m2, m3);
    float denom = Dot(m1, u);
    if (Abs(denom) < EPSILON) return 0; // Planes do not intersect in a point
    Vector3 d(p1.d, p2.d, p3.d);
    Vector3 v = Cross(m1, d);
    float ood = 1.0f / denom;
    p.x = Dot(d, u) * ood;
    p.y = Dot(m3, v) * ood;
    p.z = -Dot(m2, v) * ood;
    return 1;
}

namespace Section5 {
    // Compute the point p at which the three planes p1, p2, and p3 intersect (if at all)
    int IntersectPlanes(Plane p1, Plane p2, Plane p3, Vector3 &p)
    {
        Vector3 u = Cross(p2.n, p3.n);
        float denom = Dot(p1.n, u);
        if (Abs(denom) < EPSILON) return 0; // Planes do not intersect in a point
        p = (p1.d * u + Cross(p1.n, p3.d * p2.n - p2.d * p3.n)) / denom;
        return 1;
    }
}

// Intersect sphere s0 moving in direction d over time interval t0 <= t <= t1, against
// a stationary sphere s1. If found intersecting, return time t of collision
int TestMovingSphereSphere(Sphere s0, Vector3 d, float t0, float t1, Sphere s1, float &t)
{
    // Compute sphere bounding motion of s0 during time interval from t0 to t1
    Sphere b;
    float mid = (t0 + t1) * 0.5f;
    b.c = s0.c + d * mid;
    b.r = (mid - t0) * d.len() + s0.r;
    // If bounding sphere not overlapping s1, then no collision in this interval
    if (!TestSphereSphere(b, s1)) return 0;

    // Cannot rule collision out: recurse for more accurate testing. To terminate the
    // recursion, collision is assumed when time interval becomes sufficiently small
    if (t1 - t0 < EPSILON) {
        t = t0;
        return 1;
    }

    // Recursively test first half of interval; return collision if detected
    if (TestMovingSphereSphere(s0, d, t0, mid, s1, t)) return 1;

    // Recursively test second half of interval
    return TestMovingSphereSphere(s0, d, mid, t1, s1, t);
}

/*
//function not impl
// Test collision between objects a and b moving over the time interval
// [startTime, endTime]. When colliding, time of collision is returned in hitTime
int IntervalCollision(Object a, Object b, float startTime, float endTime, float &hitTime)
{
    // Compute the maximum distance objects a and b move over the time interval
    float maxMoveA = MaximumObjectMovementOverTime(a, startTime, endTime);
    float maxMoveB = MaximumObjectMovementOverTime(b, startTime, endTime);
    float maxMoveDistSum = maxMoveA + maxMoveB;
    // Exit if distance between a and b at start larger than sum of max movements
    float minDistStart = MinimumObjectDistanceAtTime(a, b, startTime);
    if (minDistStart > maxMoveDistSum) return 0;
    // Exit if distance between a and b at end larger than sum of max movements
    float minDistEnd = MinimumObjectDistanceAtTime(a, b, endTime);
    if (minDistEnd > maxMoveDistSum) return 0;

    // Cannot rule collision out: recurse for more accurate testing. To terminate the
    // recursion, collision is assumed when time interval becomes sufficiently small
    if (endTime - startTime < EPSILON) {
        hitTime = startTime;
        return 1;
    }
    // Recursively test first half of interval; return collision if detected
    float midTime = (startTime + endTime) * 0.5f;
    if (IntervalCollision(a, b, startTime, midTime, hitTime)) return 1;
    // Recursively test second half of interval
    return IntervalCollision(a, b, midTime, endTime, hitTime);
}
*/

// Intersect sphere s with movement vector v with plane p. If intersecting
// return time t of collision and point q at which sphere hits plane
int IntersectMovingSpherePlane(Sphere s, Vector3 v, Plane p, float &t, Vector3 &q)
{
    // Compute distance of sphere center to plane
    float dist = Dot(p.n, s.c) - p.d;
    if (Abs(dist) <= s.r) {
        // The sphere is already overlapping the plane. Set time of
        // intersection to zero and q to sphere center
        t = 0.0f;
        q = s.c;
        return 1;
    } else {
        float denom = Dot(p.n, v);
        if (denom * dist >= 0.0f) {
            // No intersection as sphere moving parallel to or away from plane
            return 0;
        } else {
            // Sphere is moving towards the plane

            // Use +r in computations if sphere in front of plane, else -r
            float r = dist > 0.0f ? s.r : -s.r;
            t = (r - dist) / denom;
            q = s.c + t * v - r * p.n;
            return 1;
        }
    }
}

// Test if sphere with radius r moving from a to b intersects with plane p
int TestMovingSpherePlane(Vector3 a, Vector3 b, float r, Plane p)
{
    // Get the distance for both a and b from plane p
    float adist = Dot(a, p.n) - p.d;
    float bdist = Dot(b, p.n) - p.d;
    // Intersects if on different sides of plane (distances have different signs)
    if (adist * bdist < 0.0f) return 1;
    // Intersects if start or end position within radius from plane
    if (Abs(adist) <= r || Abs(bdist) <= r) return 1;
    // No intersection
    return 0;
}

int TestMovingSphereSphere(Sphere s0, Sphere s1, Vector3 v0, Vector3 v1, float &t)
{
    Vector3 s = s1.c - s0.c;      // Vector3 between sphere centers
    Vector3 v = v1 - v0;          // Relative motion of s1 with respect to stationary s0
    float r = s1.r + s0.r;       // Sum of sphere radii
    float c = Dot(s, s) - r * r;
    if (c < 0.0f) {
        // Spheres initially overlapping so exit directly
        t = 0.0f;
        return 1;
    }
    float a = Dot(v, v);
    if (a < EPSILON) return 0; // Spheres not moving relative each other
    float b = Dot(v, s);
    if (b >= 0.0f) return 0;   // Spheres not moving towards each other
    float d = b * b - a * c;
    if (d < 0.0f) return 0;    // No real-valued root, spheres do not intersect

    t = (-b - std::sqrt(d)) / a;
    return 1;
}

/*
namespace Section5 {
    int TestMovingSphereSphere(Sphere s0, Sphere s1, Vector3 v0, Vector3 v1, float &t)
    {
        // Expand sphere s1 by the radius of s0
        s1.r += s0.r;
        // Subtract movement of s1 from both s0 and s1, making s1 stationary
        Vector3 v = v0 - v1;
        // Can now test directed segment s = s0.c + tv, v = (v0-v1)/||v0-v1|| against
        // the expanded sphere for intersection
        Vector3 q;
        float vlen = v.len();
        if (IntersectRaySphere(s0.c, v / vlen, s1, t, q)) {
            return t <= vlen;
        }
        return 0;
    }
}

int IntersectMovingSphereAABB(Sphere s, Vector3 d, AABB b, float &t)
{
    // Compute the AABB resulting from expanding b by sphere radius r
    AABB e = b;
    e.min.x -= s.r; e.min.y -= s.r; e.min.z -= s.r;
    e.max.x += s.r; e.max.y += s.r; e.max.z += s.r;

    // Intersect ray against expanded AABB e. Exit with no intersection if ray
    // misses e, else get intersection point p and time t as result
    Vector3 p;
    if (!IntersectRayAABB(s.c, d, e, t, p) || t > 1.0f)
        return 0;

    // Compute which min and max faces of b the intersection point p lies
    // outside of. Note, u and v cannot have the same bits set and
    // they must have at least one bit set amongst them
    int u = 0, v = 0;
    if (p.x < b.min.x) u |= 1;
    if (p.x > b.max.x) v |= 1;
    if (p.y < b.min.y) u |= 2;
    if (p.y > b.max.y) v |= 2;
    if (p.z < b.min.z) u |= 4;
    if (p.z > b.max.z) v |= 4;

    // 慜r- all set bits together into a bit mask (note: here u + v == u | v)
    int m = u + v;

    // Define line segment [c, c+d] specified by the sphere movement
    Segment seg(s.c, s.c + d);

    // If all 3 bits set (m == 7) then p is in a vertex region
    if (m == 7) {
        // Must now intersect segment [c, c+d] against the capsules of the three
        // edges meeting at the vertex and return the best time, if one or more hit
        float tmin = POS_INFINITY;
        if (IntersectSegmentCapsule(seg, Corner(b, v), Corner(b, v ^ 1), s.r, &t))
            tmin = Min(t, tmin);
        if (IntersectSegmentCapsule(seg, Corner(b, v), Corner(b, v ^ 2), s.r, &t))
            tmin = Min(t, tmin);
        if (IntersectSegmentCapsule(seg, Corner(b, v), Corner(b, v ^ 4), s.r, &t))
            tmin = Min(t, tmin);
        if (tmin == POS_INFINITY) return 0; // No intersection
        t = tmin;
        return 1; // Intersection at time t == tmin
    }
    // If only one bit set in m, then p is in a face region
    if ((m & (m - 1)) == 0) {
        // Do nothing. Time t from intersection with
        // expanded box is correct intersection time
        return 1;
    }
    // p is in an edge region. Intersect against the capsule at the edge
    return IntersectSegmentCapsule(seg, Corner(b, u ^ 7), Corner(b, v), s.r, &t);
}
*/
// Support function that returns the AABB vertex with index n
Vector3 Corner(AABB b, int n)
{
    Vector3 p;
    p.x = ((n & 1) ? b.max.x : b.min.x);
    p.y = ((n & 2) ? b.max.y : b.min.y);
    p.z = ((n & 4) ? b.max.z : b.min.z);
    return p;
}

// Intersect AABBs 慳- and 慴- moving with constant velocities va and vb.
// On intersection, return time of first and last contact in tfirst and tlast
int IntersectMovingAABBAABB(AABB a, AABB b, Vector3 va, Vector3 vb, float &tfirst, float &tlast)
{
    // Exit early if 慳- and 慴- initially overlapping
    if (TestAABBAABB(a, b)) {
        tfirst = tlast = 0.0f;
        return 1;
    }

    // Use relative velocity; effectively treating 'a' as stationary
    Vector3 v = vb - va;

    // Initialize times of first and last contact
    tfirst = 0.0f;
    tlast = 1.0f;

    // For each axis, determine times of first and last contact, if any
    for (int i = 0; i < 3; i++) {
        if (v[i] < 0.0f) {
            if (b.max[i] < a.min[i]) return 0; // Nonintersecting and moving apart
            if (a.max[i] < b.min[i]) tfirst = Max((a.max[i] - b.min[i]) / v[i], tfirst);
            if (b.max[i] > a.min[i]) tlast  = Min((a.min[i] - b.max[i]) / v[i], tlast);
        }
        if (v[i] > 0.0f) {
            if (b.min[i] > a.max[i]) return 0; // Nonintersecting and moving apart
            if (b.max[i] < a.min[i]) tfirst = Max((a.min[i] - b.max[i]) / v[i], tfirst);
            if (a.max[i] > b.min[i]) tlast = Min((a.max[i] - b.min[i]) / v[i], tlast);
        }

        // No overlap possible if time of first contact occurs after time of last contact
        if (tfirst > tlast) return 0;
    }

    return 1;
}


#endif // SECTION5_3_H
