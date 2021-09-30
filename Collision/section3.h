#ifndef SECTION3_H
#define SECTION3_H

#include "vector2.h"
#include "vector3.h"
#include "math.h"
#include <cmath>
#include <algorithm>




//----------------------------------------
//-------------------- section 3.4 --------------------
//----------------------------------------




namespace Section3 {
    //计算点p在三角形(点a,b,c)的质心坐标
    // Compute barycentric coordinates (u, v, w) for
    // point p with respect to triangle (a, b, c)
    void Barycentric(Vector3 a, Vector3 b, Vector3 c, Vector3 p, float &u, float &v, float &w)
    {
        Vector3 v0 = b - a, v1 = c - a, v2 = p - a;
        float d00 = v0.dot(v0);
        float d01 = v0.dot(v1);
        float d11 = v1.dot(v1);
        float d20 = v2.dot(v0);
        float d21 = v2.dot(v1);
        float denom = d00 * d11 - d01 * d01;
        v = (d11 * d20 - d01 * d21) / denom;
        w = (d00 * d21 - d01 * d20) / denom;
        u = 1.0f - v - w;
    }
}

inline float TriArea2D(float x1, float y1, float x2, float y2, float x3, float y3)
{
    return (x1-x2)*(y2-y3) - (x2-x3)*(y1-y2);
}

//质心坐标投影不变性
// Compute barycentric coordinates (u, v, w) for
// point p with respect to triangle (a, b, c)
void Barycentric(Vector3 a, Vector3 b, Vector3 c, Vector3 p, float &u, float &v, float &w)
{
    // Unnormalized triangle normal
    Vector3 m = (b - a).cross(c - a);
    // Nominators and one-over-denominator for u and v ratios
    float nu, nv, ood;
    // Absolute components for determining projection plane
    float x = std::abs(m.x), y = std::abs(m.y), z = std::abs(m.z);

    // Compute areas in plane of largest projection
    if (x >= y && x >= z) {
        // x is largest, project to the yz plane
        nu = TriArea2D(p.y, p.z, b.y, b.z, c.y, c.z); // Area of PBC in yz plane
        nv = TriArea2D(p.y, p.z, c.y, c.z, a.y, a.z); // Area of PCA in yz plane
        ood = 1.0f / m.x;                             // 1/(2*area of ABC in yz plane)
    } else if (y >= x && y >= z) {
        // y is largest, project to the xz plane
        nu = TriArea2D(p.x, p.z, b.x, b.z, c.x, c.z);
        nv = TriArea2D(p.x, p.z, c.x, c.z, a.x, a.z);
        ood = 1.0f / -m.y;
    } else {
        // z is largest, project to the xy plane
        nu = TriArea2D(p.x, p.y, b.x, b.y, c.x, c.y);
        nv = TriArea2D(p.x, p.y, c.x, c.y, a.x, a.y);
        ood = 1.0f / m.z;
    }
    u = nu * ood;
    v = nv * ood;
    w = 1.0f - u - v;
}

// Test if point p is contained in triangle (a, b, c)
int TestPointTriangle(Vector3 p, Vector3 a, Vector3 b, Vector3 c)
{
    float u, v, w;
    Barycentric(a, b, c, p, u, v, w);
    return v >= 0.0f && w >= 0.0f && (v + w) <= 1.0f;
}



//----------------------------------------
//-------------------- section 3.6 --------------------
//----------------------------------------



struct Plane {
    Vector3 n;  // Plane normal. Points x on the plane satisfy Dot(n,x) = d
    float d;   // d = dot(n,p) for a given point p on the plane
};

// Given three noncollinear points (ordered ccw), compute the plane equation
Plane ComputePlane(Vector3 a, Vector3 b, Vector3 c)
{
    Plane p;
    p.n = ((b - a).cross(c - a));
    p.n.normalize();
    p.d = p.n.dot(a);
    return p;
}



//----------------------------------------
//-------------------- section 3.7 --------------------
//----------------------------------------



// Test if quadrilateral (a, b, c, d) is convex
int IsConvexQuad(Vector3 a, Vector3 b, Vector3 c, Vector3 d)
{
    // Quad is nonconvex if Dot(Cross(bd, ba), Cross(bd, bc)) >= 0
    Vector3 bda = (d - b).cross(a - b);
    Vector3 bdc = (d - b).cross(c - b);
    if (bda.dot(bdc) >= 0.0f) return 0;
    // Quad is now convex iff Dot(Cross(ac, ad), Cross(ac, ab)) < 0
    Vector3 acd = (c - a).cross(d - a);
    Vector3 acb = (c - a).cross(b - a);
    return acd.dot(acb) < 0.0f;
}



//----------------------------------------
//-------------------- section 3.9 --------------------
//----------------------------------------



//给定点a,b,点集p中距离边ab最远的一点
// Return index i of point p[i] farthest from the edge ab, to the left of the edge
int PointFarthestFromEdge(Vector2 a, Vector2 b, Vector2 p[], int n)
{
    // Create edge vector and vector (counterclockwise) perpendicular to it
    Vector2 e = b - a, eperp = Vector2(-e.y, e.x);

    // Track index, 慸istance?and 憆ightmostness?of currently best point
    int bestIndex = -1;
    float maxVal = -POS_INFINITY, rightMostVal = -POS_INFINITY;

    // Test all points to find the one farthest from edge ab on the left side
    for (int i = 1; i < n; i++) {
        float d = (p[i] - a).dot(eperp); // d is proportional to distance along eperp
        float r = (p[i] - a).dot(e);     // r is proportional to distance along e
        if (d > maxVal || (d == maxVal && r > rightMostVal)) {
            bestIndex = i;
            maxVal = d;
            rightMostVal = r;
        }
    }
    return bestIndex;
}

#endif // 3_H
