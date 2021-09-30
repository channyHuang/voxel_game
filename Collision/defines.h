#ifndef DEFINES_H
#define DEFINES_H

#include "vector3.h"

typedef int int32;

struct Object {
    Object *pNextObject; // Embedded link to next hgrid object
    Vector3 pos;           // x, y (and z) position for sphere (or top left AABB corner)
    float radius;        // Radius for bounding sphere (or width of AABB)
    int bucket;          // Index of hash bucket object is in
    int level;           // Grid level for the object
    //...                  // Object data

    Vector3 center;        // Center point for object
};

#endif // DEFINES_H
