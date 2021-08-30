#pragma once

#include <vector>
#include <functional>
#include <algorithm>
#include <array>
#include <numeric>
#include <unordered_map>

#include "vector3i.h"

class TriMesh {
public:
    std::vector<Vector3> vertices;
    std::vector<Vector3> normals;
    std::vector<int> indices;

    void add_face(const Vector3i &index) {
        indices.emplace_back(index.x);
        indices.emplace_back(index.y);
        indices.emplace_back(index.z);
    }
};

class SurfaceNets {
public:
    SurfaceNets() {}
    ~SurfaceNets() {}

    TriMesh surfaceNets(std::function<float(float x, float y, float z)> const& sdfFunction,
                        Vector3i vsize,
                    float const isovalue = 0.f);

private:
    std::size_t const edges[12][2] =
    {
        { 0u, 1u },
        { 1u, 2u },
        { 2u, 3u },
        { 3u, 0u },
        { 4u, 5u },
        { 5u, 6u },
        { 6u, 7u },
        { 7u, 4u },
        { 0u, 4u },
        { 1u, 5u },
        { 2u, 6u },
        { 3u, 7u }
    };

    int const quad_neighbors[3][3] =
    {
        { 0, 1, 2 },
        { 0, 5, 4 },
        { 2, 3, 4 }
    };

    std::array<int, 3> const quad_neighbor_orders[2] =
    {
        { 0, 1, 2 },
        { 2, 1, 0 }
    };
};
