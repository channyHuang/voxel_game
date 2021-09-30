#ifndef SECTION89_H
#define SECTION89_H

#include "vector3.h"
#include "section3.h"
#include "section5.h"
typedef unsigned short uint16
typedef unsigned int uint32;


//----------------------------------------
//-------------------- section 8 --------------------
//----------------------------------------



// Constructs BSP tree from an input vector of polygons. Pass 'depth' as 0 on entry
BSPNode *BuildBSPTree(std::vector<Polygon *> &polygons, int depth)
{
    // Return NULL tree if there are no polygons
    if (polygons.empty()) return NULL;

    // Get number of polygons in the input vector
    int numPolygons = polygons.size();

    // If criterion for a leaf is matched, create a leaf node from remaining polygons
    if (depth >= MAX_DEPTH || numPolygons <= MIN_LEAF_SIZE) || ...etc...)
        return new BSPNode(polygons);

    // Select best possible partitioning plane based on the input geometry
    Plane splitPlane = PickSplittingPlane(polygons);

    std::vector<Polygon *> frontList, backList;

    // Test each polygon against the dividing plane, adding them
    // to the front list, back list, or both, as appropriate
    for (int i = 0; i < numPolygons; i++) {
        Polygon *poly = polygons[i], *frontPart, *backPart;
        switch (ClassifyPolygonToPlane(poly, splitPlane)) {
        case COPLANAR_WITH_PLANE:
            // What's done in this case depends on what type of tree is being
            // built. For a node-storing tree, the polygon is stored inside
            // the node at this level (along with all other polygons coplanar
            // with the plane). Here, for a leaf-storing tree, coplanar polygons
            // are sent to either side of the plane. In this case, to the front
            // side, by falling through to the next case
        case IN_FRONT_OF_PLANE:
            frontList.push_back(poly);
            break;
        case BEHIND_PLANE:
            backList.push_back(poly);
            break;
        case STRADDLING_PLANE:
            // Split polygon to plane and send a part to each side of the plane
            SplitPolygon(*poly, splitPlane, &frontPart, &backPart);
            frontList.push_back(frontPart);
            backList.push_back(backPart);
            break;
        }
    }

    // Recursively build child subtrees and return new tree root combining them
    BSPNode *frontTree = BuildBSPTree(frontList, depth + 1);
    BSPNode *backTree = BuildBSPTree(backList, depth + 1);
    return new BSPNode(frontTree, backTree);
}



// Given a vector of polygons, attempts to compute a good splitting plane
Plane PickSplittingPlane(std::vector<Polygon *> &polygons)
{
    // Blend factor for optimizing for balance or splits (should be tweaked)
    const float K = 0.8f;
    // Variables for tracking best splitting plane seen so far
    Plane bestPlane;
    float bestScore = FLT_MAX;

    // Try the plane of each polygon as a dividing plane
    for (int i = 0; i < polygons.size(); i++) {
        int numInFront = 0, numBehind = 0, numStraddling = 0;
        Plane plane = GetPlaneFromPolygon(polygons[i]);
        // Test against all other polygons
        for (int j = 0; j < polygons.size(); j++) {
            // Ignore testing against self
            if (i == j) continue;
            // Keep standing count of the various poly-plane relationships
            switch (ClassifyPolygonToPlane(polygons[j], plane)) {
            case POLYGON_COPLANAR_WITH_PLANE:
                /* Coplanar polygons treated as being in front of plane */
            case POLYGON_IN_FRONT_OF_PLANE:
                numInFront++;
                break;
            case POLYGON_BEHIND_PLANE:
                numBehind++;
                break;
            case POLYGON_STRADDLING_PLANE:
                numStraddling++;
                break;
            }
        }
        // Compute score as a weighted combination (based on K, with K in range
        // 0..1) between balance and splits (lower score is better)
        float score = K * numStraddling + (1.0f - K) * Abs(numInFront - numBehind);
        if (score < bestScore) {
            bestScore = score;
            bestPlane = plane;
        }
    }
    return bestPlane;
}



// Classify point p to a plane thickened by a given thickness epsilon
int ClassifyPointToPlane(Point p, Plane plane) {
    // Compute signed distance of point from plane
    float dist = Dot(plane.n, p) - plane.d;
    // Classify p based on the signed distance
    if (dist > PLANE_THICKNESS_EPSILON)
        return POINT_IN_FRONT_OF_PLANE;
    if (dist < -PLANE_THICKNESS_EPSILON)
        return POINT_BEHIND_PLANE;
    return POINT_ON_PLANE;
}


// Return value specifying whether the polygon 'poly' lies in front of,
// behind of, on, or straddles the plane 'plane'
int ClassifyPolygonToPlane(Polygon *poly, Plane plane)
{
    // Loop over all polygon vertices and count how many vertices
    // lie in front of and how many lie behind of the thickened plane
    int numInFront = 0, numBehind = 0;
    int numVerts = poly->NumVertices();
    for (int i = 0; i < numVerts; i++) {
        Point p = poly->GetVertex(i);
        switch (ClassifyPointToPlane(p, plane)) {
        case POINT_IN_FRONT_OF_PLANE:
            numInFront++;
            break;
        case POINT_BEHIND_PLANE:
            numBehind++;
            break;
        }
    }
    // If vertices on both sides of the plane, the polygon is straddling
    if (numBehind != 0 && numInFront != 0)
        return POLYGON_STRADDLING_PLANE;
    // If one or more vertices in front of the plane and no vertices behind
    // the plane, the polygon lies in front of the plane
    if (numInFront != 0)
        return POLYGON_IN_FRONT_OF_PLANE;
    // Ditto, the polygon lies behind the plane if no vertices in front of
    // the plane, and one or more vertices behind the plane
    if (numBehind != 0)
        return POLYGON_BEHIND_PLANE;
    // All vertices lie on the plane so the polygon is coplanar with the plane
    return POLYGON_COPLANAR_WITH_PLANE;
}


void SplitPolygon(Polygon &poly, Plane plane, Polygon **frontPoly, Polygon **backPoly) {
    int numFront = 0, numBack = 0;
    Point frontVerts[MAX_POINTS], backVerts[MAX_POINTS];

    // Test all edges (a, b) starting with edge from last to first vertex
    int numVerts = poly.NumVertices();
    Point a = poly.GetVertex(numVerts ?1);
    int aSide = ClassifyPointToPlane(a, plane);

    // Loop over all edges given by vertex pair (n-1, n)
    for (int n = 0; n < numVerts; n++) {
        Point b = poly.GetVertex(n);
        int bSide = ClassifyPointToPlane(b, plane);
        if (bSide == POINT_IN_FRONT_OF_PLANE) {
            if (aSide == POINT_BEHIND_PLANE) {
                // Edge (a, b) straddles, output intersection point to both sides
                Point i = IntersectEdgeAgainstPlane(a, b, plane);
                assert(ClassifyPointToPlane(i, plane) == POINT_ON_PLANE);
                frontVerts[numFront++] = backVerts[numBack++] = i;
            }
            // In all three cases, output b to the front side
            frontVerts[numFront++] = b;
        } else if (bSide == POINT_BEHIND_PLANE) {
            if (aSide == POINT_IN_FRONT_OF_PLANE) {
                // Edge (a, b) straddles plane, output intersection point
                Point i = IntersectEdgeAgainstPlane(a, b, plane);
                assert(ClassifyPointToPlane(i, plane) == POINT_ON_PLANE);
                frontVerts[numFront++] = backVerts[numBack++] = i;
            } else if (aSide == POINT_ON_PLANE) {
                // Output a when edge (a, b) goes from 憃n?to 慴ehind?plane
                backVerts[numBack++] = a;
            }
            // In all three cases, output b to the back side
            backVerts[numBack++] = b;
        } else {
            // b is on the plane. In all three cases output b to the front side
            frontVerts[numFront++] = b;
            // In one case, also output b to back side
            if (aSide == POINT_BEHIND_PLANE)
                backVerts[numBack++] = b;
        }
        // Keep b as the starting point of the next edge
        a = b;
        aSide = bSide;
    }

    // Create (and return) two new polygons from the two vertex lists
    *frontPoly = new Polygon(numFront, frontVerts);
    *backPoly = new Polygon(numBack, backVerts);
}

// Render node-storing BSP tree back-to-front w/ respect to cameraPos
void RenderBSP(BSPNode *node, Point cameraPos)
{
    // Get index of which child to visit first (0 = front, 1 = back)
    int index = ClassifyPointToPlane(cameraPos, node->plane) == POINT_IN_FRONT_OF_PLANE;

    // First visit the side the camera is NOT on
    if (node->child[index]) RenderBSP(node->child[index], cameraPos);
    // Render all polygons stored in the node
    DrawFrontfacingPolygons(node->pPolyList);
    // Then visit the other side (the one the camera is on)
    if (node->child[index ^ 1]) RenderBSP(node->child[index ^ 1], cameraPos);
}


int PointInSolidSpace(BSPNode *node, Point p)
{
    while (!node->IsLeaf()) {
        // Compute distance of point to dividing plane
        float dist = Dot(node->plane.n, p) - node->plane.d;
        if (dist > EPSILON) {
            // Point in front of plane, so traverse front of tree
            node = node->child[0];
        } else if (dist < -EPSILON) {
            // Point behind of plane, so traverse back of tree
            node = node->child[1];
        } else {
            // Point on dividing plane; must traverse both sides
            int front = PointInSolidSpace(node->child[0], p);
            int back = PointInSolidSpace(node->child[1], p);
            // If results agree, return that, else point is on boundary
            return (front == back) ? front : POINT_ON_BOUNDARY;
        }
    }
    // Now at a leaf, inside/outside status determined by solid flag
    return node->IsSolid() ? POINT_INSIDE : POINT_OUTSIDE;
}


int PointInSolidSpace(BSPNode *node, Point p)
{
    while (!node->IsLeaf()) {
        // Compute distance of point to dividing plane
        float dist = Dot(node->plane.n, p) - node->plane.d;
        // Traverse front of tree when point in front of plane, else back of tree
        node = node->child[dist <= EPSILON];
    }
    // Now at a leaf, inside/outside status determined by solid flag
    return node->IsSolid() ? POINT_INSIDE : POINT_OUTSIDE;
}


// Intersect ray/segment R(t) = p + t*d, tmin <= t <= tmax, against bsp tree
// 'node', returning time thit of first intersection with a solid leaf, if any
int RayIntersect(BSPNode *node, Point p, Vector d, float tmin, float tmax, float *thit)
{
    std::stack<BSPNode *> nodeStack;
    std::stack<float> timeStack;

    assert(node != NULL);

    while (1) {
        if (!node->IsLeaf()) {
            float denom = Dot(node->plane.n, d);
            float dist = node->plane.d - Dot(node->plane.n, p);
            int nearIndex = dist > 0.0f;
            // If denom is zero, ray runs parallel to plane. In this case,
            // just fall through to visit the near side (the one p lies on)
            if (denom != 0.0f) {
                float t = dist / denom;
                if (0.0f <= t && t <= tmax) {
                    if (t >= tmin) {
                        // Straddling, push far side onto stack, then visit near side
                        nodeStack.push(node->child[1 ^ nearIndex]);
                        timeStack.push(tmax);
                        tmax = t;
                    } else nearIndex = 1 ^ nearIndex; // 0 <= t < tmin, visit far side
                }
            }
            node = node->child[nearIndex];
        } else {
            // Now at a leaf. If it is solid, there抯 a hit at time tmin, so exit
            if (node->IsSolid()) {
                *thit = tmin;
                return 1;
            }
            // Exit if no more subtrees to visit, else pop off a node and continue
            if (nodeStack.empty()) break;
            tmin = tmax;
            node = nodeStack.top(); nodeStack.pop();
            tmax = timeStack.top(); timeStack.pop();
        }
    }

    // No hit
    return 0;
}

// Intersect polygon 'p' against the solid-leaf BSP tree 'node'
int PolygonInSolidSpace(Polygon *p, BSPNode *node)
{
    Polygon *frontPart, *backPart;
    while (!node->IsLeaf()) {
        switch (ClassifyPolygonToPlane(p, node->plane)) {
        case POLYGON_IN_FRONT_OF_PLANE:
            node = node->child[0];
            break;
        case POLYGON_BEHIND_PLANE:
            node = node->child[1];
            break;
        case POLYGON_STRADDLING_PLANE:
            SplitPolygon(*p, node->plane, &frontPart, &backPart);
            if (PolygonInSolidSpace(frontPart, node->child[0])) return 1;
            if (PolygonInSolidSpace(backPart, node->child[1])) return 1;
            // No collision
            return 0;
        }
    }
    // Now at a leaf, inside/outside status determined by solid flag
    return node->IsSolid();
}



//----------------------------------------
//-------------------- section 9 --------------------
//----------------------------------------



// Given a set s of vertices, compute a maximal set of independent vertices
Set IndependentSet(Set s)
{
    // Initialize i to the empty set
    Set i = EmptySet();
    // Loop over all vertices in the input set
    for (all vertices v in s) {
        // If unmarked and has 8 or fewer neighboring vertices...
        if (!Marked(v) && Degree(v) <= 8) {
            // Add v to the independent set and mark all of v's neighbors
            i.Add(v);
            s.MarkAllVerticesAdjacentToVertex(v);
        }
    }
    return i;
}



//----------------------------------------
//-------------------- section 11 --------------------
//----------------------------------------




// Test if segment AB intersects plane p. If so, return 1, along with
// the intersection t value and the intersection point Q. If not, return 0
int IntersectSegmentPlane(Vector3 a, Vector3 b, Plane p, float &t, Vector3 &q)
{
    // Compute t value at which the directed line ab intersects the plane
    Vector3 ab = b - a;
    t = (p.d - Dot(p.n, a)) / Dot(p.n, ab);

    // If t in [0..1] compute and return intersection point
    if (t >= 0.0f && t <= 1.0f) {
        q = a + t * ab;
        return 1;
    }
    // Else t is +INF, -INF, NaN, or not in [0..1], so no intersection
    return 0;
}

typedef uint32 uint64[2];

void Uadd64(uint64 x, uint64 y, uint64 res)
{
    uint32 a = x[1] + y[1]; // Compute sum of higher 32 bits
    uint32 b = x[0] + y[0]; // Compute sum of lower 32 bits
    if (b < x[0]) a++;      // Carry if low sum overflowed
    res[0] = b;
    res[1] = a;
}

void Umult32to64(uint32 x, uint32 y, uint64 res)
{
    uint16 xh = x >> 16, xl = x & 0xffff;
    uint16 yh = y >> 16, yl = y & 0xffff;
    uint32 a = xh * yh;
    uint32 b = xh * yl;
    uint32 c = xl * yh;
    uint32 d = xl * yl;
    d = d + (b << 16);
    if (d < (b << 16)) a++;
    d = d + (c << 16);
    if (d < (c << 16)) a++;
    a = a + (b >> 16) + (c >> 16);
    res[0] = d;
    res[1] = a;
}

// Compare rational numbers a/b and c/d
int Order(int a, int b, int c, int d)
{
    // Make c and d be nonnegative
    if (c < 0) b = -b, c = -c;
    if (d < 0) a = -a, d = -d;

    // Handle a and/or b being negative
    if (a < 0 && b < 0) {
        int olda = a, oldb = b;
        a = c; b = d; c = -olda; d = -oldb;
    }
    if (a < 0) return LESS_THAN;
    if (b < 0) return GREATER_THAN;

    // Make a <= b, exit if order becomes known
    if (a > b) {
        if (c < d) return GREATER_THAN;
        int olda = a, oldb = b;
        a = d; b = c; c = oldb; d = olda;
    }
    if (c > d) return LESS_THAN;

    // Do continued fraction expansion (given that 0<=a<=b, 0<=c<=d)
    while (a != 0 && c != 0) {
        int m = d / c;
        int n = b / a;
        if (m != n) {
            if (m < n) return LESS_THAN;
            if (m > n) return GREATER_THAN;
        }
        int olda = a, oldb = b;
        a = d % c; b = c; c = oldb % olda; d = olda;
    }
    if (a == 0) return c == 0 ? EQUAL : LESS_THAN;
    return GREATER_THAN;
}

// Test if segment ab intersects plane p. If so, return 1, along with
// an adjusted intersection point q. If not, return 0
int TestSegmentPlane(Vector3 a, Vector3 b, Plane p, Vector3 &q)
{
    // Compute t value, t=tnom/tdenom, for directed segment ab intersecting plane p
    Vector ab = b - a;
    int64 tnom = p.d - Dot(p.n, a);
    int64 tdenom = Dot(p.n, ab);

    // Exit if segment is parallel to plane
    if (tdenom == 0) return 0;

    // Ensure denominator is positive so it can be multiplied through throughout
    if (tdenom < 0) {
        tnom = -tnom;
        tdenom = -tdenom;
    }

    // If t not in [0..1], no intersection
    if (tnom < 0 || tnom > tdenom) return 0;

    // Line segment is definitely intersecting plane. Compute vector d to adjust
    // the computation of q, biasing the result to lie on the side of point a
    Vector d(0,0,0);
    int64 k = tdenom ?1;
    // If a lies behind plane p, round division other way
    if (tdenom > 0) k = -k;
    if (p.n.x > 0) d.x = k; else if (p.n.x < 0) d.x = -k;
    if (p.n.y > 0) d.y = k; else if (p.n.y < 0) d.y = -k;
    if (p.n.z > 0) d.z = k; else if (p.n.z < 0) d.z = -k;

    // Compute and return adjusted intersection point
    q = a + (tnom * ab + d) / tdenom;
    return 1;
}


#endif // SECTION89_H
