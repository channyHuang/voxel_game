#ifndef SECTION13_H
#define SECTION13_H

struct X {                  struct Y {                       struct Z {
    int8 a;                     int8 a, pad_a[7];                int64 b;
    int64 b;                    int64 b;                         int64 e;
    int8 c;                     int8 c, pad_c[1];                float f;
    int16 d;                    int16 d, pad_d[2];               int16 d;
    int64 e;                    int64 e;                         int8 a;
    float f;                    float f, pad_f[1];               int8 c;
};                          };                               };



struct S {
    int32 value;
    int32 count;
    ...
} elem[1000];

// Find the index of the element with largest value
int index = 0;
for (int i = 0; i < 1000; i++)
    if (elem[i].value > elem[index].value) index = i;
// Increment the count field of that element
elem[index].count++;



// Hot fields of S
struct S1 {
    int32 value;
} elem1[1000];

// Cold fields of S
struct S2 {
    int32 count;
    ...
} elem2[1000];

// Find the index of the element with largest value
int index = 0;
for (int i = 0; i < 1000; i++)
    if (elem1[i].value > elem1[index].value) index = i;
// Increment the count field of that element
elem2[index].count++;


// Loop through and process all 4n elements
for (int i = 0; i < 4 * n; i++)
    Process(elem[i]);



const int kLookAhead = 4;            // Experiment to find which fixed distance works best
Prefetch(&elem[0]);                  // Prefetch first four elements
for (int i = 0; i < 4 * n; i += 4) { // Unrolled. Process and step 4 elements at a time
    Prefetch(&elem[i + kLookAhead]); // Prefetch cache line a few elements ahead
    Process(elem[i + 0]);            // Process the elements that have already
    Process(elem[i + 1]);            // been fetched into memory on the
    Process(elem[i + 2]);            // previous iteration
    Process(elem[i + 3]);
}



void PreorderTraversal(Node *pNode)
{
    Prefetch(pNode->left);           // Greedily prefetch left traversal path
    Process(pNode);                  // Process the current node
    Prefetch(pNode->right);          // Greedily prefetch right traversal path
    PreorderTraversal(pNode->left);  // Recursively visit left subtree
    PreorderTraversal(pNode->right); // then recursively visit right subtree
}



Elem a = elem[0];          // Load the cache line of the first four array elements
for (int i = 0; i < 4 * n; i += 4) {
    Elem e = elem[i + 4];  // Cache miss. Fetches next cache line using nonblocking load
    Elem b = elem[i + 1];  // These following three values are in cache and...
    Elem c = elem[i + 2];  // ...are loaded as hit-under-miss; no stalls
    Elem d = elem[i + 3];
    Process(a);            // Process the data from the cache line that...
    Process(b);            // ...was fetched in the previous iteration
    Process(c);
    Process(d);
    a = e;                 // e now in cache, and so is b, c, and d of next iteration
}

union KDNode {
    float splitVal_type;  // nonleaf, type 00 = x, 01 = y, 10 = z-split
    int32 leafIndex_type; // leaf, type 11
};



// Align tree root to start of cache line (64-byte aligned)
void ComputeChildPointers(KDNode *pRoot, KDNode *pNode, KDNode **pLeft, KDNode **pRight)
{
    int32 nodeAddress = (int32)pNode;
    int32 nodeIndex = (nodeAddress & 0x3f) >> 2; // node index within cache line (0-14)
    int32 leftAddress, rightAddress;
    if (nodeIndex < 7) {
        // Three out of four, children are at 2n+1 and 2n+2 within current cache line
        leftAddress = nodeAddress + (nodeAddress & 0x3f) + 4;
        rightAddress = leftAddress + 4;
    } else {
        // The children are roots of subtrees at some other cache lines
        int32 rootAddress = (int32)pRoot;
        int32 cacheLineFromStart = (nodeAddress - rootAddress) >> 6;
        int32 leftCacheLine = 16 * cacheLineFromStart + 1;
        int32 bitIndex = nodeIndex - 7; // (0-7)
        leftAddress = rootAddress + leftCacheLine * 64 + bitIndex * 2 * 64;
        rightAddress = leftAddress + 64; // at next consecutive cache line
    }
    *pLeft = (KDNode *)leftAddress;
    *pRight = (KDNode *)rightAddress;
}



// Align tree root to start of cache line (64-byte aligned)
void ComputeChildPointers(KDNode *pRoot, KDNode *pNode, KDNode **pLeft, KDNode **pRight)
{
    int32 nodeAddress = (int32)pNode;
    int32 nodeIndex = (nodeAddress & 0x3f) >> 2; // node index within cache line (0-14)
    int32 leftAddress, rightAddress;
    if (nodeIndex < 7) {
        // Three out of four children are at 2n+1 and 2n+2 within current cache line
        leftAddress = nodeAddress + (nodeAddress & 0x3f) + 4;
        rightAddress = leftAddress + 4;
    } else {
        // The children are roots of subtrees at some other cache lines
        int32 rootAddress = (int32)pRoot;
        int32 bitIndex = nodeIndex - 7; // (0-7)
        // Last word on cache line specifies linking of subtrees
        int32 linkWord = *((int32 *)(nodeAddress | 0x3c));
        assert(linkWord & (1 << bitIndex)); // must be set
        int32 offset = PopCount8(linkWord & ((1 << bitIndex) - 1));
        leftAddress = rootAddress + ((linkWord >> 8) + offset * 2) * 64;
        rightAddress = leftAddress + 64; // at next consecutive cache line
    }
    *pLeft = (KDNode *)leftAddress;
    *pRight = (KDNode *)rightAddress;
}



inline int32 PopCount8(int32 n)
{
    n = n - ((n & 0xaa) >> 1);          // Form partial sums of two bits each
    n = (n & 0x33) + ((n >> 2) & 0x33); // Form partial sums of four bits each
    return (n + (n >> 4)) & 0x0f;       // Add four-bit sums together and mask result
}


// Uncompressed AABB
struct AABB {
    float & operator [](int n) { return ((float *)this)[n]; }
    Vector3 min;          // The minimum point defining the AABB
    Vector3 max;          // The maximum point defining the AABB
};

// Compressed AABB node
struct PackedAABBNode {
    uint8 flags;        // Bits determining if new extent belong to left or right child
    uint8 newExtent[6]; // New extents, quantized within the parent volume
    ...                 // Potential links to children nodes
};

void ComputeChildAABBs(AABB &parent, PackedAABBNode &node, AABB *left, AABB *right)
{
    for (int i = 0; i < 6; i++) {
        int xyz = i & 3;
        // Compute the actual side value from the quantized newExtent[] value
        float min = parent.min[xyz], max = parent.max[xyz];
        float val = min + (max - min) * (float)node.newExtent[i] / 255.0f;
        // Test bits to see which child gets parent side value and which get new value
        if (node.flags & (1 << i)) {
            (*left)[i] = parent[i];  // Left child gets parent's bound
            (*right)[i] = val;       // ...right child gets computed bound
        } else {
            (*left)[i] = val;        // Left child gets computed bound
            (*right)[i] = parent[i]; // ...right child gets parent's bound
        }
    }
}


// Array of all vertices, to avoid repeating vertices between triangles
Vector3 vertex[NUM_VERTICES];

// Array of all triangles, to avoid repeating triangles between leaves
struct IndexedTriangle {
    int16 i0, i1, i2;
} triangle[NUM_TRIANGLES];

// Variable-size leaf node structure, containing indices into triangle array
struct Leaf {
    int16 numTris; // number of triangles in leaf node
    int16 tri[1];  // first element of variable-size array, with one or more triangles
};



// Explicitly declared triangle, no indirection
struct Triangle {
    Vector3 p0, p1, p2;
};

#define NUM_LEAFCACHE_ENTRIES (128-1)          // 2^N-1 for fast hash key masking w/ '&'

// Cached leaf node, no indirection
struct CachedLeaf {
    Leaf *pLeaf;                               // hash key corresponding to leaf id
    int16 numTris;                             // number of triangles in leaf node
    Triangle triangle[MAX_TRIANGLES_PER_LEAF]; // the cached triangles
} leafCache[NUM_LEAFCACHE_ENTRIES];



// Compute direct-mapped hash table index from pLeaf pointer (given as input)
int32 hashIndex = ((int32)pLeaf >> 2) & NUM_LEAFCACHE_ENTRIES;

// If leaf not in cache at this address, cache it
CachedLeaf *pCachedLeaf = &leafCache[hashIndex];
if (pCachedLeaf->pLeaf != pLeaf) {
    // Set cache key to be pLeaf pointer
    pCachedLeaf->pLeaf = pLeaf;
    // Fetch all triangles and store in cached leaf (linearized)
    int numTris = pLeaf->numTris;
    pCachedLeaf->numTris = numTris;
    for (int i = 0; i < numTris; i++) {
        Triangle *pCachedTri = &pCachedLeaf->triangle[i];
        IndexedTriangle *pTri = &triangle[pLeaf->tri[i]];
        pCachedTri->p0 = vertex[pTri->i0];
        pCachedTri->p1 = vertex[pTri->i1];
        pCachedTri->p2 = vertex[pTri->i2];
    }
}

// Leaf now in cache so use cached leaf node data for processing leaf node
int numTris = pCachedLeaf->numTris;
for (int i = 0; i < numTris; i++) {
    Triangle *pCachedTri = &pCachedLeaf->triangle[i];
    DoSomethingWithTriangle(pCachedTri->p0, pCachedTri->p1, pCachedTri->p2);
}
...


int n;
int *p1 = &n; // *p1 references the location of n
int *p2 = &n; // *p2 also references the location of n



void TransformPoint(Vector3 *pOut, Matrix &m, Vector3 &in)
{
    pOut->x = m[0][0] * in.x + m[0][1] * in.y + m[0][2] * in.z;
    pOut->y = m[1][0] * in.x + m[1][1] * in.y + m[1][2] * in.z;
    pOut->z = m[2][0] * in.x + m[2][1] * in.y + m[2][2] * in.z;
}



void TransformPoint(Vector3 *pOut, Matrix &m, Vector3 &in)
{
    float tmpx, tmpy, tmpz;
    tmpx = m[0][0] * in.x + m[0][1] * in.y + m[0][2] * in.z;
    tmpy = m[1][0] * in.x + m[1][1] * in.y + m[1][2] * in.z;
    tmpz = m[2][0] * in.x + m[2][1] * in.y + m[2][2] * in.z;
    pOut->x = tmpx; pOut->y = tmpy; pOut->z = tmpz;
}


void Foo(float *v, int *n)
{
    for (int i = 0; i < *n; i++)
        v[i] += 1.0f;
}



void Foo(float *v, int *n)
{
    int t = *n;
    for (int i = 0; i < t; i++)
        v[i] += 1.0f;
}



uint32 i;
float f;
i = *((int32 *)&f); // Illegal way of getting the integer representation of the float f



union {
    uint32 i;
    float f;
} u;
u.f = f;
uint32 i = u.i;


void TransformPoint(Vector3 * restrict pOut, Matrix & restrict m, Vector3 & restrict in)
{
    pOut->x = m[0][0] * in.x + m[0][1] * in.y + m[0][2] * in.z;
    pOut->y = m[1][0] * in.x + m[1][1] * in.y + m[1][2] * in.z;
    pOut->z = m[2][0] * in.x + m[2][1] * in.y + m[2][2] * in.z;
}



void ComplexMult(float *a, float *b, float *c)
{
    a[0] = b[0]*c[0] - b[1]*c[1]; // real part in a[0] (b[0] and c[0])
    a[1] = b[0]*c[1] + b[1]*c[0]; // imaginary part in a[1] (b[1] and c[1])
}


#endif // SECTION13_H
