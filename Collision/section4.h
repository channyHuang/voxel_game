#ifndef SECTION3_6_H
#define SECTION3_6_H

#include "math.h"
#include "vector2.h"
#include "vector3.h"
#include "matrix33.h"

//----------------------------------------
//-------------------- section 4.1 ~ 4.2 --------------------
//----------------------------------------



// region R = { (x, y, z) | min.x<=x<=max.x, min.y<=y<=max.y, min.z<=z<=max.z }
struct AABB {
    Vector3 min;
    Vector3 max;
};

int TestAABBAABB(AABB a, AABB b)
{
    // Exit with no intersection if separated along an axis
    if (a.max[0] < b.min[0] || a.min[0] > b.max[0]) return 0;
    if (a.max[1] < b.min[1] || a.min[1] > b.max[1]) return 0;
    if (a.max[2] < b.min[2] || a.min[2] > b.max[2]) return 0;
    // Overlapping on all axes means AABBs are intersecting
    return 1;
}

namespace Section4_1 {
    // region R = {(x, y, z) | min.x<=x<=min.x+dx, min.y<=y<=min.y+dy, min.z<=z<=min.z+dz }
    struct AABB {
        Vector3 min;
        float d[3];   // diameter or width extents (dx, dy, dz)
    };

    int TestAABBAABB(AABB a, AABB b)
    {
        float t;
        if ((t = a.min[0] - b.min[0]) > b.d[0] || -t > a.d[0]) return 0;
        if ((t = a.min[1] - b.min[1]) > b.d[1] || -t > a.d[1]) return 0;
        if ((t = a.min[2] - b.min[2]) > b.d[2] || -t > a.d[2]) return 0;
        return 1;
    }
}

namespace Section4_2 {
    // region R = { (x, y, z) | |c.x-x|<=rx, |c.y-y|<=ry, |c.z-z|<=rz }
    struct AABB {
        Vector3 c; // center Vector3 of AABB
        float r[3]; // radius or halfwidth extents (rx, ry, rz)
    };

    int TestAABBAABB(AABB a, AABB b)
    {
        if (std::abs(a.c[0] - b.c[0]) > (a.r[0] + b.r[0])) return 0;
        if (std::abs(a.c[1] - b.c[1]) > (a.r[1] + b.r[1])) return 0;
        if (std::abs(a.c[2] - b.c[2]) > (a.r[2] + b.r[2])) return 0;
        return 1;
    }

    namespace Section4_2_1 {
        int TestAABBAABB(AABB a, AABB b)
        {
            int r;
            r = a.r[0] + b.r[0]; if ((unsigned int)(a.c[0] - b.c[0] + r) > r + r) return 0;
            r = a.r[1] + b.r[1]; if ((unsigned int)(a.c[1] - b.c[1] + r) > r + r) return 0;
            r = a.r[2] + b.r[2]; if ((unsigned int)(a.c[2] - b.c[2] + r) > r + r) return 0;
            return 1;
        }
    }

    // Transform AABB a by the matrix m and translation t,
    // find maximum extents, and store result into AABB b.
    void UpdateAABB(AABB a, float m[3][3], float t[3], AABB &b)
    {
        for (int i = 0; i < 3; i++) {
            b.c[i] = t[i];
            b.r[i] = 0.0f;
            for (int j = 0; j < 3; j++) {
                b.c[i] += m[i][j] * a.c[j];
                b.r[i] += std::abs(m[i][j]) * a.r[j];
            }
        }
    }
}

// Returns indices imin and imax into pt[] array of the least and
// most, respectively, distant Vector3s along the direction dir
void ExtremeVector3sAlongDirection(Vector3 dir, Vector3 pt[], int n, int *imin, int *imax)
{
    float minproj = POS_INFINITY, maxproj = -POS_INFINITY;
    for (int i = 0; i < n; i++) {
        // Project vector from origin to Vector3 onto direction vector
        float proj = pt[i].dot(dir);
        // Keep track of least distant Vector3 along direction vector
        if (proj < minproj) {
            minproj = proj;
            *imin = i;
        }
        // Keep track of most distant Vector3 along direction vector
        if (proj > maxproj) {
            maxproj = proj;
            *imax = i;
        }
    }
}

// Transform AABB a by the matrix m and translation t,
// find maximum extents, and store result into AABB b.
void UpdateAABB(AABB a, float m[3][3], float t[3], AABB &b)
{
    // For all three axes
    for (int i = 0; i < 3; i++) {
        // Start by adding in translation
        b.min[i] = b.max[i] = t[i];
        // Form extent by summing smaller and larger terms respectively
        for (int j = 0; j < 3; j++) {
            float e = m[i][j] * a.min[j];
            float f = m[i][j] * a.max[j];
            if (e < f) {
                b.min[i] += e;
                b.max[i] += f;
            } else {
                b.min[i] += f;
                b.max[i] += e;
            }
        }
    }
}



//----------------------------------------
//-------------------- section 4.3 --------------------
//----------------------------------------



// Region R = { (x, y, z) | (x-c.x)^2 + (y-c.y)^2 + (z-c.z)^2 <= r^2 }
struct Sphere {
    Vector3 c; // Sphere center
    float r; // Sphere radius

    Sphere() {}
    Sphere(const Vector3 &p) { c = p; r = 0; }
    Sphere(const Vector3 &p, const Vector3 &q) {
        c = (p + q) * 0.5f;
        r = p.distanceTo(q);
    }
};

int TestSphereSphere(Sphere a, Sphere b)
{
    // Calculate squared distance between centers
    Vector3 d = a.c - b.c;
    float dist2 = d.dot(d);
    // Spheres intersect if squared distance is less than squared sum of radii
    float radiusSum = a.r + b.r;
    return dist2 <= radiusSum * radiusSum;
}

//计算包围球
// Compute indices to the two most separated Vector3s of the (up to) six Vector3s
// defining the AABB encompassing the Vector3 set. Return these as min and max.
void MostSeparatedVector3sOnAABB(int &min, int &max, Vector3 pt[], int numPts)
{
    // First find most extreme Vector3s along principal axes
    int minx = 0, maxx = 0, miny = 0, maxy = 0, minz = 0, maxz = 0;
    for (int i = 1; i < numPts; i++) {
        if (pt[i].x < pt[minx].x) minx = i;
        if (pt[i].x > pt[maxx].x) maxx = i;
        if (pt[i].y < pt[miny].y) miny = i;
        if (pt[i].y > pt[maxy].y) maxy = i;
        if (pt[i].z < pt[minz].z) minz = i;
        if (pt[i].z > pt[maxz].z) maxz = i;
    }

    // Compute the squared distances for the three pairs of Vector3s
    float dist2x = (pt[maxx] - pt[minx]).dot(pt[maxx] - pt[minx]);
    float dist2y = (pt[maxy] - pt[miny]).dot(pt[maxy] - pt[miny]);
    float dist2z = (pt[maxz] - pt[minz]).dot(pt[maxz] - pt[minz]);
    // Pick the pair (min,max) of Vector3s most distant
    min = minx;
    max = maxx;
    if (dist2y > dist2x && dist2y > dist2z) {
        max = maxy;
        min = miny;
    }
    if (dist2z > dist2x && dist2z > dist2y) {
        max = maxz;
        min = minz;
    }
}

void SphereFromDistantVector3s(Sphere &s, Vector3 pt[], int numPts)
{
    // Find the most separated Vector3 pair defining the encompassing AABB
    int min, max;
    MostSeparatedVector3sOnAABB(min, max, pt, numPts);

    // Set up sphere to just encompass these two Vector3s
    s.c = (pt[min] + pt[max]) * 0.5f;
    s.r = (pt[max] - s.c).dot(pt[max] - s.c);
    s.r = std::sqrt(s.r);
}

// Given Sphere s and Vector3 p, update s (if needed) to just encompass p
void SphereOfSphereAndPt(Sphere &s, Vector3 &p)
{
    // Compute squared distance between Vector3 and sphere center
    Vector3 d = p - s.c;
    float dist2 = d.dot(d);
    // Only update s if Vector3 p is outside it
    if (dist2 > s.r * s.r) {
        float dist = std::sqrt(dist2);
        float newRadius = (s.r + dist) * 0.5f;
        float k = (newRadius - s.r) / dist;
        s.r = newRadius;
        s.c += d * k;
    }
}

void RitterSphere(Sphere &s, Vector3 pt[], int numPts)
{
    // Get sphere encompassing two approximately most distant Vector3s
    SphereFromDistantVector3s(s, pt, numPts);

    // Grow sphere to include all Vector3s
    for (int i = 0; i < numPts; i++)
        SphereOfSphereAndPt(s, pt[i]);
}

//计算方差
// Compute variance of a set of 1D values
float Variance(float x[], int n)
{
    float u = 0.0f;
    for (int i = 0; i < n; i++)
        u += x[i];
    u /= n;
    float s2 = 0.0f;
    for (int i = 0; i < n; i++)
        s2 += (x[i] - u) * (x[i] - u);
    return s2 / n;
}

//计算协方差矩阵
void CovarianceMatrix(Matrix33 &cov, Vector3 pt[], int numPts)
{
    float oon = 1.0f / (float)numPts;
    Vector3 c = Vector3(0.0f, 0.0f, 0.0f);
    float e00, e11, e22, e01, e02, e12;

    // Compute the center of mass (centroid) of the Vector3s
    for (int i = 0; i < numPts; i++)
        c += pt[i];
    c *= oon;

    // Compute covariance elements
    e00 = e11 = e22 = e01 = e02 = e12 = 0.0f;
    for (int i = 0; i < numPts; i++) {
        // Translate Vector3s so center of mass is at origin
        Vector3 p = pt[i] - c;
        // Compute covariance of translated Vector3s
        e00 += p.x * p.x;
        e11 += p.y * p.y;
        e22 += p.z * p.z;
        e01 += p.x * p.y;
        e02 += p.x * p.z;
        e12 += p.y * p.z;
    }
    // Fill in the covariance matrix elements
    cov[0][0] = e00 * oon;
    cov[1][1] = e11 * oon;
    cov[2][2] = e22 * oon;
    cov[0][1] = cov[1][0] = e01 * oon;
    cov[0][2] = cov[2][0] = e02 * oon;
    cov[1][2] = cov[2][1] = e12 * oon;
}

//Jacobi分解，A=J^T*A*J
// 2-by-2 Symmetric Schur decomposition. Given an n-by-n symmetric matrix
// and indicies p, q such that 1 <= p < q <= n, computes a sine-cosine pair
// (s, c) that will serve to form a Jacobi rotation matrix.
//
// See Golub, Van Loan, Matrix Computations, 3rd ed, p428
void SymSchur2(Matrix33 &a, int p, int q, float &c, float &s)
{
    if (std::abs(a[p][q]) > 0.0001f) {
        float r = (a[q][q] - a[p][p]) / (2.0f * a[p][q]);
        float t;
        if (r >= 0.0f)
            t = 1.0f / (r + std::sqrt(1.0f + r*r));
        else
            t = -1.0f / (-r + std::sqrt(1.0f + r*r));
        c = 1.0f / std::sqrt(1.0f + t*t);
        s = t * c;
    } else {
        c = 1.0f;
        s = 0.0f;
    }
}

// Computes the eigenvectors and eigenvalues of the symmetric matrix A using
// the classic Jacobi method of iteratively updating A as A = J^T * A * J,
// where J = J(p, q, theta) is the Jacobi rotation matrix.
//
// On exit, v will contain the eigenvectors, and the diagonal elements
// of a are the corresponding eigenvalues.
//
// See Golub, Van Loan, Matrix Computations, 3rd ed, p428
void Jacobi(Matrix33 &a, Matrix33 &v)
{
    int i, j, n, p, q;
    float prevoff, c, s;
    Matrix33 J, b, t;

    // Initialize v to identity matrix
    for (i = 0; i < 3; i++) {
        v[i][0] = v[i][1] = v[i][2] = 0.0f;
        v[i][i] = 1.0f;
    }

    // Repeat for some maximum number of iterations
    const int MAX_ITERATIONS = 50;
    for (n = 0; n < MAX_ITERATIONS; n++) {
        // Find largest off-diagonal std::absolute element a[p][q]
        p = 0; q = 1;
        for (i = 0; i < 3; i++) {
            for (j = 0; j < 3; j++) {
                if (i == j) continue;
                if (std::abs(a[i][j]) > std::abs(a[p][q])) {
                    p = i;
                    q = j;
                }
            }
        }

        // Compute the Jacobi rotation matrix J(p, q, theta)
        // (This code can be optimized for the three different cases of rotation)
        SymSchur2(a, p, q, c, s);
        for (i = 0; i < 3; i++) {
            J[i][0] = J[i][1] = J[i][2] = 0.0f;
            J[i][i] = 1.0f;
        }
        J[p][p] =  c; J[p][q] = s;
        J[q][p] = -s; J[q][q] = c;

        // Cumulate rotations into what will contain the eigenvectors
        v = v * J;

        // Make 'a' more diagonal, until just eigenvalues remain on diagonal
        a = (J.Transpose() * a) * J;

        // Compute "norm" of off-diagonal elements
        float off = 0.0f;
        for (i = 0; i < 3; i++) {
            for (j = 0; j < 3; j++) {
                if (i == j) continue;
                off += a[i][j] * a[i][j];
            }
        }
        /* off = sqrt(off); not needed for norm comparison */

        // Stop when norm no longer decreasing
        if (n > 2 && off >= prevoff)
            return;

        prevoff = off;
    }
}

void EigenSphere(Sphere &eigSphere, Vector3 pt[], int numPts)
{
    Matrix33 m, v;

    // Compute the covariance matrix m
    CovarianceMatrix(m, pt, numPts);
    // Decompose it into eigenvectors (in v) and eigenvalues (in m)
    Jacobi(m, v);

    // Find the component with largest magnitude eigenvalue (largest spread)
    Vector3 e;
    int maxc = 0;
    float maxf, maxe = std::abs(m[0][0]);
    if ((maxf = std::abs(m[1][1])) > maxe) maxc = 1, maxe = maxf;
    if ((maxf = std::abs(m[2][2])) > maxe) maxc = 2, maxe = maxf;
    e[0] = v[0][maxc];
    e[1] = v[1][maxc];
    e[2] = v[2][maxc];

    // Find the most extreme Vector3s along direction 'e'
    int imin, imax;
    ExtremeVector3sAlongDirection(e, pt, numPts, &imin, &imax);
    Vector3 minpt = pt[imin];
    Vector3 maxpt = pt[imax];

    float dist = std::sqrt((maxpt - minpt).dot(maxpt - minpt));
    eigSphere.r = dist * 0.5f;
    eigSphere.c = (minpt + maxpt) * 0.5f;
}

void RitterEigenSphere(Sphere &s, Vector3 pt[], int numPts)
{
    // Start with sphere from maximum spread
    EigenSphere(s, pt, numPts);

    // Grow sphere to include all Vector3s
    for (int i = 0; i < numPts; i++)
        SphereOfSphereAndPt(s, pt[i]);
}

void RitterIterative(Sphere &s, Vector3 pt[], int numPts)
{
    const int NUM_ITER = 8;
    RitterSphere(s, pt, numPts);
    Sphere s2 = s;
    for (int k = 0; k < NUM_ITER; k++) {
        // Shrink sphere somewhat to make it an underestimate (not bound)
        s2.r = s2.r * 0.95f;

        // Make sphere bound data again
        for (int i = 0; i < numPts; i++) {
            // Swap pt[i] with pt[j], where j randomly from interval [i+1,numPts-1]
            //DoRandomSwap();
            SphereOfSphereAndPt(s2, pt[i]);
        }

        // Update s whenever a tighter sphere is found
        if (s2.r < s.r) s = s2;
    }
}

bool Vector3InsideSphere(const Vector3 &p, Sphere &s) {
    if (p.distanceTo(s.c) <= s.r) return true;
    return false;
}

//最小圆覆盖Welzl算法
Sphere WelzlSphere(Vector3 pt[], unsigned int numPts, Vector3 sos[], unsigned int numSos)
{
    // if no input Vector3s, the recursion has bottomed out. Now compute an
    // exact sphere based on Vector3s in set of support (zero through four Vector3s)
    if (numPts == 0) {
        switch (numSos) {
        case 0: return Sphere();
        case 1: return Sphere(sos[0]);
        case 2: return Sphere(sos[0], sos[1]);
        //case 3: return Sphere(sos[0], sos[1], sos[2]);
        //case 4: return Sphere(sos[0], sos[1], sos[2], sos[3]);
        }
    }
    // Pick a Vector3 at "random" (here just the last Vector3 of the input set)
    int index = numPts - 1;
    // Recursively compute the smallest bounding sphere of the remaining Vector3s
    Sphere smallestSphere = WelzlSphere(pt, numPts - 1, sos, numSos); // (*)
    // If the selected Vector3 lies inside this sphere, it is indeed the smallest
    if(Vector3InsideSphere(pt[index], smallestSphere))
        return smallestSphere;
    // Otherwise, update set of support to additionally contain the new Vector3
    sos[numSos] = pt[index];
    // Recursively compute the smallest sphere of remaining Vector3s with new s.o.s.
    return WelzlSphere(pt, numPts - 1, sos, numSos + 1);
}



//----------------------------------------
//-------------------- section 4.4 --------------------
//----------------------------------------



// Region R = { x | x = c+r*u[0]+s*u[1]+t*u[2] }, |r|<=e[0], |s|<=e[1], |t|<=e[2]
struct OBB {
    Vector3 c;     // OBB center Vector3
    Vector3 u[3]; // Local x-, y-, and z-axes
    Vector3 e;    // Positive halfwidth extents of OBB along each axis
};

int TestOBBOBB(OBB &a, OBB &b)
{
    float ra, rb;
    Matrix33 R, AbsR;

    // Compute rotation matrix expressing b in a's coordinate frame
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
            R[i][j] = (a.u[i]).dot(b.u[j]);

    // Compute translation vector t
    Vector3 t = b.c - a.c;
    // Bring translation into a's coordinate frame
    t = Vector3(t.dot(a.u[0]), t.dot(a.u[1]), t.dot(a.u[2]));

    // Compute common subexpressions. Add in an epsilon term to
    // counteract arithmetic errors when two edges are parallel and
    // their cross product is (near) null (see text for details)
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
            AbsR[i][j] = std::abs(R[i][j]) + EPSILON;

    // Test axes L = A0, L = A1, L = A2
    for (int i = 0; i < 3; i++) {
        ra = a.e[i];
        rb = b.e[0] * AbsR[i][0] + b.e[1] * AbsR[i][1] + b.e[2] * AbsR[i][2];
        if (std::abs(t[i]) > ra + rb) return 0;
    }

    // Test axes L = B0, L = B1, L = B2
    for (int i = 0; i < 3; i++) {
        ra = a.e[0] * AbsR[0][i] + a.e[1] * AbsR[1][i] + a.e[2] * AbsR[2][i];
        rb = b.e[i];
        if (std::abs(t[0] * R[0][i] + t[1] * R[1][i] + t[2] * R[2][i]) > ra + rb) return 0;
    }

    // Test axis L = A0 x B0
    ra = a.e[1] * AbsR[2][0] + a.e[2] * AbsR[1][0];
    rb = b.e[1] * AbsR[0][2] + b.e[2] * AbsR[0][1];
    if (std::abs(t[2] * R[1][0] - t[1] * R[2][0]) > ra + rb) return 0;

    // Test axis L = A0 x B1
    ra = a.e[1] * AbsR[2][1] + a.e[2] * AbsR[1][1];
    rb = b.e[0] * AbsR[0][2] + b.e[2] * AbsR[0][0];
    if (std::abs(t[2] * R[1][1] - t[1] * R[2][1]) > ra + rb) return 0;

    // Test axis L = A0 x B2
    ra = a.e[1] * AbsR[2][2] + a.e[2] * AbsR[1][2];
    rb = b.e[0] * AbsR[0][1] + b.e[1] * AbsR[0][0];
    if (std::abs(t[2] * R[1][2] - t[1] * R[2][2]) > ra + rb) return 0;

    // Test axis L = A1 x B0
    ra = a.e[0] * AbsR[2][0] + a.e[2] * AbsR[0][0];
    rb = b.e[1] * AbsR[1][2] + b.e[2] * AbsR[1][1];
    if (std::abs(t[0] * R[2][0] - t[2] * R[0][0]) > ra + rb) return 0;

    // Test axis L = A1 x B1
    ra = a.e[0] * AbsR[2][1] + a.e[2] * AbsR[0][1];
    rb = b.e[0] * AbsR[1][2] + b.e[2] * AbsR[1][0];
    if (std::abs(t[0] * R[2][1] - t[2] * R[0][1]) > ra + rb) return 0;

    // Test axis L = A1 x B2
    ra = a.e[0] * AbsR[2][2] + a.e[2] * AbsR[0][2];
    rb = b.e[0] * AbsR[1][1] + b.e[1] * AbsR[1][0];
    if (std::abs(t[0] * R[2][2] - t[2] * R[0][2]) > ra + rb) return 0;

    // Test axis L = A2 x B0
    ra = a.e[0] * AbsR[1][0] + a.e[1] * AbsR[0][0];
    rb = b.e[1] * AbsR[2][2] + b.e[2] * AbsR[2][1];
    if (std::abs(t[1] * R[0][0] - t[0] * R[1][0]) > ra + rb) return 0;

    // Test axis L = A2 x B1
    ra = a.e[0] * AbsR[1][1] + a.e[1] * AbsR[0][1];
    rb = b.e[0] * AbsR[2][2] + b.e[2] * AbsR[2][0];
    if (std::abs(t[1] * R[0][1] - t[0] * R[1][1]) > ra + rb) return 0;

    // Test axis L = A2 x B2
    ra = a.e[0] * AbsR[1][2] + a.e[1] * AbsR[0][2];
    rb = b.e[0] * AbsR[2][1] + b.e[1] * AbsR[2][0];
    if (std::abs(t[1] * R[0][2] - t[0] * R[1][2]) > ra + rb) return 0;

    // Since no separating axis found, the OBBs must be intersecting
    return 1;
}

// Compute the center point, 'c', and axis orientation, u[0] and u[1], of
// the minimum area rectangle in the xy plane containing the points pt[].
float MinAreaRect(Vector2 pt[], int numPts, Vector2 &c, Vector2 u[2])
{
    float minArea = POS_INFINITY;

    // Loop through all edges; j trails i by 1, modulo numPts
    for (int i = 0, j = numPts - 1; i < numPts; j = i, i++) {
        // Get current edge e0 (e0x,e0y), normalized
        Vector2 e0 = pt[i] - pt[j];
        e0 /= e0.len();

        // Get an axis e1 orthogonal to edge e0
        Vector2 e1 = Vector2(-e0.y, e0.x); // = Perp2D(e0)

        // Loop through all points to get maximum extents
        float min0 = 0.0f, min1 = 0.0f, max0 = 0.0f, max1 = 0.0f;
        for (int k = 0; k < numPts; k++) {
            // Project points onto axes e0 and e1 and keep track
            // of minimum and maximum values along both axes
            Vector2 d = pt[k] - pt[j];
            float dot = d.dot(e0);
            if (dot < min0) min0 = dot;
            if (dot > max0) max0 = dot;
            dot = d.dot(e1);
            if (dot < min1) min1 = dot;
            if (dot > max1) max1 = dot;
        }
        float area = (max0 - min0) * (max1 - min1);

        // If best so far, remember area, center, and axes
        if (area < minArea) {
            minArea = area;
            c = pt[j] + 0.5f * ((min0 + max0) * e0 + (min1 + max1) * e1);
            u[0] = e0; u[1] = e1;
        }
    }
    return minArea;
}


//----------------------------------------
//-------------------- section 4.5 --------------------
//----------------------------------------



// Region R = { x | (x - [a + (b - a)*t])^2 <= r }, 0 <= t <= 1
struct Capsule {
    Vector3 a;      // Medial line segment start point
    Vector3 b;      // Medial line segment end point
    float r;      // Radius
};

// Region R = { x | (x - [a + u[0]*s + u[1]*t])^2 <= r }, 0 <= s,t <= 1
struct Lozenge {
    Vector3 a;      // Origin
    Vector3 u[2];  // The two edges axes of the rectangle
    float r;      // Radius
};

//SqDistPointSegment and ClosestPtSegmentSegment are in section5
int TestSphereCapsule(Sphere s, Capsule capsule)
{
    // Compute (squared) distance between sphere center and capsule line segment
    float dist2 = 0;//SqDistPointSegment(capsule.a, capsule.b, s.c);

    // If (squared) distance smaller than (squared) sum of radii, they collide
    float radius = s.r + capsule.r;
    return dist2 <= radius * radius;
}

int TestCapsuleCapsule(Capsule capsule1, Capsule capsule2)
{
    // Compute (squared) distance between the inner structures of the capsules
    float s, t;
    Vector3 c1, c2;
    float dist2 = 0;//ClosestPtSegmentSegment(capsule1.a, capsule1.b, capsule2.a, capsule2.b, s, t, c1, c2);

    // If (squared) distance smaller than (squared) sum of radii, they collide
    float radius = capsule1.r + capsule2.r;
    return dist2 <= radius * radius;
}



//----------------------------------------
//-------------------- section 4.6 --------------------
//----------------------------------------



// Region R = { (x, y, z) | dNear <= a*x + b*y + c*z <= dFar }
struct Slab {
    float n[3];  // Normal n = (a, b, c)
    float dNear; // Signed distance from origin for near plane (dNear)
    float dFar;  // Signed distance from origin for far plane (dFar)
};

struct DOP8 {
    float min[4]; // Minimum distance (from origin) along axes 0 to 3
    float max[4]; // Maximum distance (from origin) along axes 0 to 3
};

template<int k> struct KDOP {
    float min[k];
    float max[k];
};

/*
int TestKDOPKDOP(KDOP &a, KDOP &b, int k)
{
    // Check if any intervals are non-overlapping, return if so
    for (int i = 0; i < k / 2; i++)
        if (a.min[i] > b.max[i] || a.max[i] < b.min[i])
            return 0;

    // All intervals are overlapping, so k-DOPs must intersect
    return 1;
}
*/

// Compute 8-DOP for object vertices v[] in world space
// using the axes (1,1,1), (1,1,-1), (1,-1,1) and (-1,1,1)
void ComputeDOP8(Vector3 v[], int numPts, DOP8 &dop8)
{
    // Initialize 8-DOP to an empty volume
    dop8.min[0] = dop8.min[1] = dop8.min[2] = dop8.min[3] = POS_INFINITY;
    dop8.max[0] = dop8.max[1] = dop8.max[2] = dop8.max[3] = -POS_INFINITY;

    // For each point, update 8-DOP bounds if necessary
    float value;
    for (int i = 0; i < numPts; i++) {
        // Axis 0 = (1,1,1)
        value = v[i].x + v[i].y + v[i].z;
        if (value < dop8.min[0]) dop8.min[0] = value;
        else if (value > dop8.max[0]) dop8.max[0] = value;

        // Axis 1 = (1,1,-1)
        value = v[i].x + v[i].y - v[i].z;
        if (value < dop8.min[1]) dop8.min[1] = value;
        else if (value > dop8.max[1]) dop8.max[1] = value;

        // Axis 2 = (1,-1,1)
        value = v[i].x - v[i].y + v[i].z;
        if (value < dop8.min[2]) dop8.min[2] = value;
        else if (value > dop8.max[2]) dop8.max[2] = value;

        // Axis 3 = (-1,1,1)
        value = -v[i].x + v[i].y + v[i].z;
        if (value < dop8.min[3]) dop8.min[3] = value;
        else if (value > dop8.max[3]) dop8.max[3] = value;
    }
}


#endif // SECTION3_6_H
