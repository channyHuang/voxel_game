#ifndef SECTION5_H
#define SECTION5_H

#include <vector>
#include "section4.h"

//self add
float Dot(const Vector3 &a, const Vector3 &b) {
    return a.dot(b);
}
Vector3 Cross(const Vector3 &a, const Vector3 &b) {
    return a.cross(b);
}
float Abs(const float v) {
    return std::abs(v);
}
Vector3 Normalize(const Vector3 &a) {
    return a.normalizedCopy();
}
float Max(const float a, const float b) {
    return std::max(a, b);
}
float Max(const float a, const float b, const float c) {
    return std::max(std::max(a, b), c);
}
float Min(const float a, const float b) {
    return std::min(a, b);
}
float Min(const float a, const float b, const float c) {
    return std::min(std::min(a, b), c);
}
float Swap(float &a, float &b) {
    float t = a;
    a = b;
    b = t;
}
float Sqr(float a) {
    return a * a;
}


//section3
struct Plane {
    Vector3 n;  // Plane normal. Vector3s x on the plane satisfy Dot(n,x) = d
    float d;   // d = dot(n,p) for a given point p on the plane
};



//----------------------------------------
//-------------------- section 5.1 --------------------
//----------------------------------------



//点q距平面p最近的点
Vector3 ClosestPtPointPlane(Vector3 q, Plane p)
{
    float t = (p.n.dot(q) - p.d) / p.n.dot(p.n);
    return q - t * p.n;
}

namespace Section5 {
    //单位法向量
    Vector3 ClosestPtPointPlane(Vector3 q, Plane p)
    {
        float t = p.n.dot(q) - p.d;
        return q - t * p.n;
    }
}

//点q距平面p的有向距离
float DistPointPlane(Vector3 q, Plane p)
{
    // return Dot(q, p.n) - p.d; if plane equation normalized (||p.n||==1)
    return (p.n.dot(q) - p.d) / p.n.dot(p.n);
}

// Given segment ab and point c, computes closest point d on ab.
// Also returns t for the position of d, d(t) = a + t*(b - a)
void ClosestPtPointSegment(Vector3 c, Vector3 a, Vector3 b, float &t, Vector3 &d)
{
    Vector3 ab = b - a;
    // Project c onto ab, computing parameterized position d(t) = a + t*(b ?a)
    t = (c - a).dot(ab) / ab.dot(ab);
    // If outside segment, clamp t (and therefore d) to the closest endpoint
    if (t < 0.0f) t = 0.0f;
    if (t > 1.0f) t = 1.0f;
    // Compute projected position from the clamped t
    d = a + t * ab;
}

namespace {
    // Given segment ab and point c, computes closest point d on ab.
    // Also returns t for the parametric position of d, d(t) = a + t*(b - a)
    void ClosestPtPointSegment(Vector3 c, Vector3 a, Vector3 b, float &t, Vector3 &d)
    {
        Vector3 ab = b - a;
        // Project c onto ab, but deferring divide by Dot(ab, ab)
        t = (c - a).dot(ab);
        if (t <= 0.0f) {
            // c projects outside the [a,b] interval, on the a side; clamp to a
            t = 0.0f;
            d = a;
        } else {
            float denom = ab.dot(ab); // Always nonnegative since denom = ||ab||^2
            if (t >= denom) {
                // c projects outside the [a,b] interval, on the b side; clamp to b
                t = 1.0f;
                d = b;
            } else {
                // c projects inside the [a,b] interval; must do deferred divide now
                t = t / denom;
                d = a + t * ab;
            }
        }
    }
}

// Returns the squared distance between point c and segment ab
float SqDistPointSegment(Vector3 a, Vector3 b, Vector3 c)
{
    Vector3 ab = b - a, ac = c - a, bc = c - b;
    float e = ac.dot(ab);
    // Handle cases where c projects outside ab
    if (e <= 0.0f) return ac.dot(ac);
    float f = ab.dot(ab);
    if (e >= f) return bc.dot(bc);
    // Handle case where c projects onto ab
    return ac.dot(ac) - e * e / f;
}

// Given point p, return the point q on or in AABB b, that is closest to p
void ClosestPtPointAABB(Vector3 p, AABB b, Vector3 &q)
{
    // For each coordinate axis, if the point coordinate value is
    // outside box, clamp it to the box, else keep it as is
    for (int i = 0; i < 3; i++) {
        float v = p[i];
        if (v < b.min[i]) v = b.min[i]; // v = max(v, b.min[i])
        if (v > b.max[i]) v = b.max[i]; // v = min(v, b.max[i])
        q[i] = v;
    }
}

// Computes the square distance between a point p and an AABB b
float SqDistPointAABB(Vector3 p, AABB b)
{
    float sqDist = 0.0f;
    for (int i = 0; i < 3; i++) {
        // For each axis count any excess distance outside box extents
        float v = p[i];
        if (v < b.min[i]) sqDist += (b.min[i] - v) * (b.min[i] - v);
        if (v > b.max[i]) sqDist += (v - b.max[i]) * (v - b.max[i]);
    }
    return sqDist;
}

// Given point p, return point q on (or in) OBB b, closest to p
void ClosestPtPointOBB(Vector3 p, OBB b, Vector3 &q)
{
    Vector3 d = p - b.c;
    // Start result at center of box; make steps from there
    q = b.c;
    // For each OBB axis...
    for (int i = 0; i < 3; i++) {
        // ...project d onto that axis to get the distance
        // along the axis of d from the box center
        float dist = d.dot(b.u[i]);
        // If distance farther than the box extents, clamp to the box
        if (dist > b.e[i]) dist = b.e[i];
        if (dist < -b.e[i]) dist = -b.e[i];
        // Step that distance along the axis to get world coordinate
        q += dist * b.u[i];
    }
}

// Computes the square distance between point p and OBB b
float SqDistPointOBB(Vector3 p, OBB b)
{
    Vector3 closest;
    ClosestPtPointOBB(p, b, closest);
    float sqDist = (closest - p).dot(closest - p);
    return sqDist;
}

namespace {
    // Computes the square distance between point p and OBB b
    float SqDistPointOBB(Vector3 p, OBB b)
    {
        Vector3 v = p - b.c;
        float sqDist = 0.0f;
        for (int i = 0; i < 3; i++) {
            // Project vector from box center to p on each axis, getting the distance
            // of p along that axis, and count any excess distance outside box extents
            float d = v.dot(b.u[i]), excess = 0.0f;
            if (d < -b.e[i])
                excess = d + b.e[i];
            else if (d > b.e[i])
                excess = d - b.e[i];
            sqDist += excess * excess;
        }
        return sqDist;
    }
}

struct Rect {
    Vector3 c;     // center point of rectangle
    Vector3 u[2]; // unit vectors determining local x and y axes for the rectangle
    float e[2];  // the halfwidth extents of the rectangle along the axes
};

// Given point p, return point q on (or in) Rect r, closest to p
void ClosestPtPointRect(Vector3 p, Rect r, Vector3 &q)
{
    Vector3 d = p - r.c;
    // Start result at center of rect; make steps from there
    q = r.c;
    // For each rect axis...
    for (int i = 0; i < 2; i++) {
        // ...project d onto that axis to get the distance
        // along the axis of d from the rect center
        float dist = d.dot(r.u[i]);
        // If distance farther than the rect extents, clamp to the rect
        if (dist > r.e[i]) dist = r.e[i];
        if (dist < -r.e[i]) dist = -r.e[i];
        // Step that distance along the axis to get world coordinate
        q += dist * r.u[i];
    }
}

namespace {
    // Return point q on (or in) rect (specified by a, b, and c), closest to given point p
    void ClosestPtPointRect(Vector3 p, Vector3 a, Vector3 b, Vector3 c, Vector3 &q)
    {
        Vector3 ab = b - a; // vector across rect
        Vector3 ac = c - a; // vector down rect
        Vector3 d = p - a;
        // Start result at top-left corner of rect; make steps from there
        q = a;
        // Clamp p?(projection of p to plane of r) to rectangle in the across direction
        float dist = d.dot(ab);
        float maxdist = ab.dot(ab);
        if (dist >= maxdist)
            q += ab;
        else if (dist > 0.0f)
            q += (dist / maxdist) * ab;
        // Clamp p?(projection of p to plane of r) to rectangle in the down direction
        dist = d.dot(ac);
        maxdist = ac.dot(ac);
        if (dist >= maxdist)
            q += ac;
        else if (dist > 0.0f)
            q += (dist / maxdist) * ac;
    }
}

Vector3 ClosestPtPointTriangle(Vector3 p, Vector3 a, Vector3 b, Vector3 c)
{
    Vector3 ab = b - a;
    Vector3 ac = c - a;
    Vector3 bc = c - b;

    // Compute parametric position s for projection P' of P on AB,
    // P' = A + s*AB, s = snom/(snom+sdenom)
    float snom = (p - a).dot(ab), sdenom = (p - b).dot(a - b);

    // Compute parametric position t for projection P' of P on AC,
    // P' = A + t*AC, s = tnom/(tnom+tdenom)
    float tnom = (p - a).dot(ac), tdenom = (p - c).dot(a - c);

    if (snom <= 0.0f && tnom <= 0.0f) return a; // Vertex region early out

    // Compute parametric position u for projection P' of P on BC,
    // P' = B + u*BC, u = unom/(unom+udenom)
    float unom = (p - b).dot(bc), udenom = (p - c).dot(b - c);

    if (sdenom <= 0.0f && unom <= 0.0f) return b; // Vertex region early out
    if (tdenom <= 0.0f && udenom <= 0.0f) return c; // Vertex region early out


    // P is outside (or on) AB if the triple scalar product [N PA PB] <= 0
    Vector3 n = (b - a).cross(c - a);
    float vc = n.dot((a - p).cross(b - p));
    // If P outside AB and within feature region of AB,
    // return projection of P onto AB
    if (vc <= 0.0f && snom >= 0.0f && sdenom >= 0.0f)
        return a + snom / (snom + sdenom) * ab;

    // P is outside (or on) BC if the triple scalar product [N PB PC] <= 0
    float va = n.dot((b - p).cross(c - p));
    // If P outside BC and within feature region of BC,
    // return projection of P onto BC
    if (va <= 0.0f && unom >= 0.0f && udenom >= 0.0f)
        return b + unom / (unom + udenom) * bc;

    // P is outside (or on) CA if the triple scalar product [N PC PA] <= 0
    float vb = n.dot((c - p).cross(a - p));
    // If P outside CA and within feature region of CA,
    // return projection of P onto CA
    if (vb <= 0.0f && tnom >= 0.0f && tdenom >= 0.0f)
        return a + tnom / (tnom + tdenom) * ac;

    // P must project inside face region. Compute Q using barycentric coordinates
    float u = va / (va + vb + vc);
    float v = vb / (va + vb + vc);
    float w = 1.0f - u - v; // = vc / (va + vb + vc)
    return u * a + v * b + w * c;
}

namespace Section5 {
    Vector3 ClosestPtPointTriangle(Vector3 p, Vector3 a, Vector3 b, Vector3 c)
    {
        // Check if P in vertex region outside A
        Vector3 ab = b - a;
        Vector3 ac = c - a;
        Vector3 ap = p - a;
        float d1 = ab.dot(ap);
        float d2 = ac.dot(ap);
        if (d1 <= 0.0f && d2 <= 0.0f) return a; // barycentric coordinates (1,0,0)

        // Check if P in vertex region outside B
        Vector3 bp = p - b;
        float d3 = ab.dot(bp);
        float d4 = ac.dot(bp);
        if (d3 >= 0.0f && d4 <= d3) return b; // barycentric coordinates (0,1,0)

        // Check if P in edge region of AB, if so return projection of P onto AB
        float vc = d1*d4 - d3*d2;
        if (vc <= 0.0f && d1 >= 0.0f && d3 <= 0.0f) {
            float v = d1 / (d1 - d3);
            return a + v * ab; // barycentric coordinates (1-v,v,0)
        }

        // Check if P in vertex region outside C
        Vector3 cp = p - c;
        float d5 = ab.dot(cp);
        float d6 = ac.dot(cp);
        if (d6 >= 0.0f && d5 <= d6) return c; // barycentric coordinates (0,0,1)

        // Check if P in edge region of AC, if so return projection of P onto AC
        float vb = d5*d2 - d1*d6;
        if (vb <= 0.0f && d2 >= 0.0f && d6 <= 0.0f) {
            float w = d2 / (d2 - d6);
            return a + w * ac; // barycentric coordinates (1-w,0,w)
        }

        // Check if P in edge region of BC, if so return projection of P onto BC
        float va = d3*d6 - d5*d4;
        if (va <= 0.0f && (d4 - d3) >= 0.0f && (d5 - d6) >= 0.0f) {
            float w = (d4 - d3) / ((d4 - d3) + (d5 - d6));
            return b + w * (c - b); // barycentric coordinates (0,1-w,w)
        }

        // P inside face region. Compute Q through its barycentric coordinates (u,v,w)
        float denom = 1.0f / (va + vb + vc);
        float v = vb * denom;
        float w = vc * denom;
        return a + ab * v + ac * w; // = u*a + v*b + w*c, u = va * denom = 1.0f - v - w
    }
}

// Test if point p lies outside plane through abc
int PointOutsideOfPlane(Vector3 p, Vector3 a, Vector3 b, Vector3 c)
{
    return Dot(p - a, Cross(b - a, c - a)) >= 0.0f; // [AP AB AC] >= 0
}

namespace {
    // Test if point p and d lie on opposite sides of plane through abc
    int PointOutsideOfPlane(Vector3 p, Vector3 a, Vector3 b, Vector3 c, Vector3 d)
    {
        float signp = Dot(p - a, Cross(b - a, c - a)); // [AP AB AC]
        float signd = Dot(d - a, Cross(b - a, c - a)); // [AD AB AC]
        // Points on opposite sides if expression signs are opposite
        return signp * signd < 0.0f;
    }
}

Vector3 ClosestPtPointTetrahedron(Vector3 p, Vector3 a, Vector3 b, Vector3 c, Vector3 d)
{
    // Start out assuming point inside all halfspaces, so closest to itself
    Vector3 closestPt = p;
    float bestSqDist = POS_INFINITY;
    // If point outside face abc then compute closest point on abc
    if (PointOutsideOfPlane(p, a, b, c)) {
        Vector3 q = ClosestPtPointTriangle(p, a, b, c);
        float sqDist = Dot(q - p, q - p);
        // Update best closest point if (squared) distance is less than current best
        if (sqDist < bestSqDist) bestSqDist = sqDist, closestPt = q;
    }
    // Repeat test for face acd
    if (PointOutsideOfPlane(p, a, c, d)) {
        Vector3 q = ClosestPtPointTriangle(p, a, c, d);
        float sqDist = Dot(q - p, q - p);
        if (sqDist < bestSqDist) bestSqDist = sqDist, closestPt = q;
    }
    // Repeat test for face adb
    if (PointOutsideOfPlane(p, a, d, b)) {
        Vector3 q = ClosestPtPointTriangle(p, a, d, b);
        float sqDist = Dot(q - p, q - p);
        if (sqDist < bestSqDist) bestSqDist = sqDist, closestPt = q;
    }
    // Repeat test for face bdc
    if (PointOutsideOfPlane(p, b, d, c)) {
        Vector3 q = ClosestPtPointTriangle(p, b, d, c);
        float sqDist = Dot(q - p, q - p);
        if (sqDist < bestSqDist) bestSqDist = sqDist, closestPt = q;
    }
    return closestPt;
}

// Clamp n to lie within the range [min, max]
float Clamp(float n, float min, float max) {
    if (n < min) return min;
    if (n > max) return max;
    return n;
}

// Computes closest points C1 and C2 of S1(s)=P1+s*(Q1-P1) and
// S2(t)=P2+t*(Q2-P2), returning s and t. Function result is squared
// distance between between S1(s) and S2(t)
float ClosestPtSegmentSegment(Vector3 p1, Vector3 q1, Vector3 p2, Vector3 q2,
                              float &s, float &t, Vector3 &c1, Vector3 &c2)
{
    Vector3 d1 = q1 - p1; // Direction vector of segment S1
    Vector3 d2 = q2 - p2; // Direction vector of segment S2
    Vector3 r = p1 - p2;
    float a = Dot(d1, d1); // Squared length of segment S1, always nonnegative
    float e = Dot(d2, d2); // Squared length of segment S2, always nonnegative
    float f = Dot(d2, r);

    // Check if either or both segments degenerate into points
    if (a <= EPSILON && e <= EPSILON) {
        // Both segments degenerate into points
        s = t = 0.0f;
        c1 = p1;
        c2 = p2;
        return Dot(c1 - c2, c1 - c2);
    }
    if (a <= EPSILON) {
        // First segment degenerates into a point
        s = 0.0f;
        t = f / e; // s = 0 => t = (b*s + f) / e = f / e
        t = Clamp(t, 0.0f, 1.0f);
    } else {
        float c = Dot(d1, r);
        if (e <= EPSILON) {
            // Second segment degenerates into a point
            t = 0.0f;
            s = Clamp(-c / a, 0.0f, 1.0f); // t = 0 => s = (b*t - c) / a = -c / a
        } else {
            // The general nondegenerate case starts here
            float b = Dot(d1, d2);
            float denom = a*e-b*b; // Always nonnegative

            // If segments not parallel, compute closest point on L1 to L2, and
            // clamp to segment S1. Else pick arbitrary s (here 0)
            if (denom != 0.0f) {
                s = Clamp((b*f - c*e) / denom, 0.0f, 1.0f);
            } else s = 0.0f;

#if 1
            // Compute point on L2 closest to S1(s) using
            // t = Dot((P1+D1*s)-P2,D2) / Dot(D2,D2) = (b*s + f) / e
            t = (b*s + f) / e;

            // If t in [0,1] done. Else clamp t, recompute s for the new value
            // of t using s = Dot((P2+D2*t)-P1,D1) / Dot(D1,D1)= (t*b - c) / a
            // and clamp s to [0, 1]
            if (t < 0.0f) {
                t = 0.0f;
                s = Clamp(-c / a, 0.0f, 1.0f);
            } else if (t > 1.0f) {
                t = 1.0f;
                s = Clamp((b - c) / a, 0.0f, 1.0f);
            }
#else
            float tnom = b*s + f;
            if (tnom < 0.0f) {
                t = 0.0f;
                s = Clamp(-c / a, 0.0f, 1.0f);
            } else if (tnom > e) {
                t = 1.0f;
                s = Clamp((b - c) / a, 0.0f, 1.0f);
            } else {
                t = tnom / e;
            }
#endif
        }
    }

    c1 = p1 + d1 * s;
    c2 = p2 + d2 * t;
    return Dot(c1 - c2, c1 - c2);
}

// Returns 2 times the signed triangle area. The result is positive if
// abc is ccw, negative if abc is cw, zero if abc is degenerate.
float Signed2DTriArea(Vector3 a, Vector3 b, Vector3 c)
{
    return (a.x - c.x) * (b.y - c.y) - (a.y - c.y) * (b.x - c.x);
}

// Test if segments ab and cd overlap. If they do, compute and return
// intersection t value along ab and intersection position p
int Test2DSegmentSegment(Vector3 a, Vector3 b, Vector3 c, Vector3 d, float &t, Vector3 &p)
{
    // Sign of areas correspond to which side of ab points c and d are
    float a1 = Signed2DTriArea(a, b, d); // Compute winding of abd (+ or -)
    float a2 = Signed2DTriArea(a, b, c); // To intersect, must have sign opposite of a1

    // If c and d are on different sides of ab, areas have different signs
    if (a1 * a2 < 0.0f) {
        // Compute signs for a and b with respect to segment cd
        float a3 = Signed2DTriArea(c, d, a); // Compute winding of cda (+ or -)
        // Since area is constant a1-a2 = a3-a4, or a4=a3+a2-a1
//      float a4 = Signed2DTriArea(c, d, b); // Must have opposite sign of a3
        float a4 = a3 + a2 - a1;
        // Points a and b on different sides of cd if areas have different signs
        if (a3 * a4 < 0.0f) {
            // Segments intersect. Find intersection point along L(t)=a+t*(b-a).
            // Given height h1 of a over cd and height h2 of b over cd,
            // t = h1 / (h1 - h2) = (b*h1/2) / (b*h1/2 - b*h2/2) = a3 / (a3 - a4),
            // where b (the base of the triangles cda and cdb, i.e., the length
            // of cd) cancels out.
            t = a3 / (a3 - a4);
            p = a + t * (b - a);
            return 1;
        }
    }

    // Segments not intersecting (or collinear)
    return 0;
}

#endif // SECTION5_H
