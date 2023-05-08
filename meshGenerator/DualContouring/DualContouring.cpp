#include "DualContouring.h"

namespace DualContouring {
    namespace {
        static const float TRANSITION_CELL_SCALE = 0.25;

        // Multiply integer math results by this
        static const float FIXED_FACTOR = 1.f / 256.f;

        inline float tof(int8_t v)
        {
            return static_cast<float>(v) / 256.f;
        }

        inline int8_t tos(uint8_t v)
        {
            return v - 128;
        }

        // Values considered negative have a sign bit of 1
        inline uint8_t sign(int8_t v)
        {
            return (v >> 7) & 1;
        }

        // Wrapped to invert SDF data, Transvoxel apparently works backwards?
        inline uint8_t get_voxel(const VoxelBuffer &vb, int x, int y, int z, int channel)
        {
            return 255 - static_cast<uint8_t>(vb.get_voxel(x, y, z, channel));
        }

        inline uint8_t get_voxel(const VoxelBuffer &vb, const Vector3i &pos, int channel)
        {
            return get_voxel(vb, pos.x, pos.y, pos.z, channel);
        }

        inline Vector3 normalized_not_null(const Vector3 &n)
        {
            float lengthsq = n.lenSqr();
            if (lengthsq == 0)
            {
                return Vector3(0, 1, 0);
            }
            else
            {
                float length = std::sqrt(lengthsq);
                return Vector3(n.x / length, n.y / length, n.z / length);
            }
        }


        Vector3 clamp_to(const Vector3 &v, const Vector3i &vmin, const Vector3i &vmax) {
            Vector3 res = v;
            if (v.x < vmin.x) res.x = vmin.x;
            if (v.y < vmin.y) res.y = vmin.y;
            if (v.z < vmin.z) res.z = vmin.z;
            if (v.x > vmax.x) res.x = vmax.x;
            if (v.y > vmax.y) res.y = vmax.y;
            if (v.z > vmax.z) res.z = vmax.z;
            return res;
        }

    };

    std::vector<Vector3> DualContouring::getNormal(const VoxelBuffer &voxels, Vector3i pos) {
        std::vector<Vector3> vNormals(8, Vector3(0, 0, 0));

        const Vector3i min_pos = Vector3i(0);

        std::array<int8_t, 8> cell_samples;
        std::array<Vector3, 8> corner_gradients;
        std::array<Vector3i, 8> corner_positions;
        //    6-------7
        //   /|      /|
        //  / |     / |  Corners
        // 4-------5  |
        // |  2----|--3
        // | /     | /   z y
        // |/      |/    |/
        // 0-------1     o--x
        //
        // Warning: temporarily includes padding. It is undone later.
        corner_positions[0] = Vector3i(pos.x, pos.y, pos.z);
        corner_positions[1] = Vector3i(pos.x + 1, pos.y, pos.z);
        corner_positions[2] = Vector3i(pos.x, pos.y + 1, pos.z);
        corner_positions[3] = Vector3i(pos.x + 1, pos.y + 1, pos.z);
        corner_positions[4] = Vector3i(pos.x, pos.y, pos.z + 1);
        corner_positions[5] = Vector3i(pos.x + 1, pos.y, pos.z + 1);
        corner_positions[6] = Vector3i(pos.x, pos.y + 1, pos.z + 1);
        corner_positions[7] = Vector3i(pos.x + 1, pos.y + 1, pos.z + 1);

        // Get the value of cells.
        // Negative values are "solid" and positive are "air".
        // Due to raw cells being unsigned 8-bit, they get converted to signed.
        for (unsigned int i = 0; i < corner_positions.size(); ++i)
        {
            cell_samples[i] = tos(get_voxel(voxels, corner_positions[i], VoxelBuffer::CHANNEL_SDF));
        }

        // Concatenate the sign of cell values to obtain the case code.
        // Index 0 is the less significant bit, and index 7 is the most significant bit.
        uint8_t case_code = sign(cell_samples[0]);
        case_code |= (sign(cell_samples[1]) << 1);
        case_code |= (sign(cell_samples[2]) << 2);
        case_code |= (sign(cell_samples[3]) << 3);
        case_code |= (sign(cell_samples[4]) << 4);
        case_code |= (sign(cell_samples[5]) << 5);
        case_code |= (sign(cell_samples[6]) << 6);
        case_code |= (sign(cell_samples[7]) << 7);

        if (case_code == 0 || case_code == 255)
        {
            // If the case_code is 0 or 255, there is no triangulation to do
            return vNormals;
        }

        // Compute normals
        for (unsigned int i = 0; i < corner_positions.size(); ++i)
        {
            Vector3i p = corner_positions[i];

            float nx = tof(tos(get_voxel(voxels, p.x - 1, p.y, p.z, VoxelBuffer::CHANNEL_SDF)));
            float ny = tof(tos(get_voxel(voxels, p.x, p.y - 1, p.z, VoxelBuffer::CHANNEL_SDF)));
            float nz = tof(tos(get_voxel(voxels, p.x, p.y, p.z - 1, VoxelBuffer::CHANNEL_SDF)));
            float px = tof(tos(get_voxel(voxels, p.x + 1, p.y, p.z, VoxelBuffer::CHANNEL_SDF)));
            float py = tof(tos(get_voxel(voxels, p.x, p.y + 1, p.z, VoxelBuffer::CHANNEL_SDF)));
            float pz = tof(tos(get_voxel(voxels, p.x, p.y, p.z + 1, VoxelBuffer::CHANNEL_SDF)));

            //get_gradient_normal(nx, px, ny, py, nz, pz, cell_samples[i]);
            corner_gradients[i] = Vector3(nx - px, ny - py, nz - pz);

            // Undo padding here. From this point, corner positions are actual positions.
            corner_positions[i] = (corner_positions[i] - min_pos) << 0;
        }

        // For cells occurring along the minimal boundaries of a block,
        // the preceding cells needed for vertex reuse may not exist.
        // In these cases, we allow new vertex creation on additional edges of a cell.
        // While iterating through the cells in a block, a 3-bit mask is maintained whose bits indicate
        // whether corresponding bits in a direction code are valid

        uint8_t direction_validity_mask =
            (pos.x > min_pos.x ? 1 : 0) |
            ((pos.y > min_pos.y ? 1 : 0) << 1) |
            ((pos.z > min_pos.z ? 1 : 0) << 2);

        uint8_t regular_cell_class_index = Transvoxel::get_regular_cell_class(case_code);
        const Transvoxel::RegularCellData &regular_cell_data = Transvoxel::get_regular_cell_data(regular_cell_class_index);
        uint8_t triangle_count = regular_cell_data.geometryCounts & 0x0f;
        uint8_t vertex_count = (regular_cell_data.geometryCounts & 0xf0) >> 4;

        std::array<int, 12> cell_vertex_indices;
        for (int i = 0; i < 12; ++i)
        {
            cell_vertex_indices[i] = -1;
        }

        // For each vertex in the case
        for (unsigned int i = 0; i < vertex_count; ++i)
        {
            // The case index maps to a list of 16-bit codes providing information about the edges on which the vertices lie.
            // The low byte of each 16-bit code contains the corner indexes of the edge��s endpoints in one nibble each,
            // and the high byte contains the mapping code shown in Figure 3.8(b)
            unsigned short rvd = Transvoxel::get_regular_vertex_data(case_code, i);
            uint8_t edge_code_low = rvd & 0xff;
            uint8_t edge_code_high = (rvd >> 8) & 0xff;

            // Get corner indexes in the low nibble (always ordered so the higher comes last)
            uint8_t v0 = (edge_code_low >> 4) & 0xf;
            uint8_t v1 = edge_code_low & 0xf;

            // Get voxel values at the corners
            int sample0 = cell_samples[v0]; // called d0 in the paper
            int sample1 = cell_samples[v1]; // called d1 in the paper

            // Get interpolation position
            // We use an 8-bit fraction, allowing the new vertex to be located at one of 257 possible
            // positions  along  the  edge  when  both  endpoints  are included.
            int t = (sample1 << 8) / (sample1 - sample0);

            float t0 = static_cast<float>(t) / 256.f;
            float t1 = static_cast<float>(0x100 - t) / 256.f;
            int ti0 = t;
            int ti1 = 0x100 - t;

            const Vector3i p0 = corner_positions[v0];
            const Vector3i p1 = corner_positions[v1];

            if (t & 0xff)
            {
                // Vertex is between p0 and p1 (inside the edge)

                // Each edge of a cell is assigned an 8-bit code, as shown in Figure 3.8(b),
                // that provides a mapping to a preceding cell and the coincident edge on that preceding cell
                // for which new vertex creation  was  allowed.
                // The high nibble of this code indicates which direction to go in order to reach the correct preceding cell.
                // The bit values 1, 2, and 4 in this nibble indicate that we must subtract one
                // from the x, y, and/or z coordinate, respectively.
                uint8_t reuse_dir = (edge_code_high >> 4) & 0xf;
                uint8_t reuse_vertex_index = edge_code_high & 0xf;

                // TODO Some re-use opportunities are missed on negative sides of the block,
                // but I don't really know how to fix it...
                // You can check by "shaking" every vertex randomly in a shader based on its index,
                // you will see vertices touching the -X, -Y or -Z sides of the block aren't connected

                bool present = (reuse_dir & direction_validity_mask) == reuse_dir;

                //if (!present || cell_vertex_indices[i] == -1)
                {
                    Vector3i primary = p0 * ti0 + p1 * ti1;
                    Vector3 primaryf = primary.to_vec3() * FIXED_FACTOR;
                    Vector3 normal = normalized_not_null(corner_gradients[v0] * t0 + corner_gradients[v1] * t1);

                    normal.normalize();
                    vNormals[v0] += normal;
                    vNormals[v0].normalize();
                    vNormals[v1] += normal;
                    vNormals[v1].normalize();
                }
            }
            else if (t == 0 && v1 == 7)
            {
                Vector3i primary = p1;
                Vector3 primaryf = primary.to_vec3();
                Vector3 normal = normalized_not_null(corner_gradients[v1]);

                normal.normalize();
                vNormals[v1] += normal;
                vNormals[v1].normalize();
                vNormals[v0] += normal;
                vNormals[v0].normalize();
            }
            else
            {
                uint8_t reuse_dir = (t == 0 ? v1 ^ 7 : v0 ^ 7);
                bool present = (reuse_dir & direction_validity_mask) == reuse_dir;

                //if (!present || cell_vertex_indices[i] < 0)
                {
                    Vector3i primary = t == 0 ? p1 : p0;
                    Vector3 primaryf = primary.to_vec3();
                    Vector3 normal = normalized_not_null(corner_gradients[t == 0 ? v1 : v0]);

                    normal.normalize();
                    vNormals[(t == 0 ? v1 : v0)] += normal;
                    vNormals[(t == 0 ? v1 : v0)].normalize();
                    vNormals[(t == 0 ? v0 : v1)] += normal;
                    vNormals[(t == 0 ? v0 : v1)].normalize();
                }
            }
        }
        return vNormals;
    }

    void DualContouring::constructVertices(const VoxelBuffer &buffer, TriMesh &mesh, VoxelBuffer &vertIdx,
        const int min_padding, const int blocksize_with_padding) {
        Vector3i vSize = buffer.get_size() - Vector3i(1, 1, 1);
        // Reused for speed (no dynamic allocations in inner loop):
        std::vector<Plane>  planes;
        std::vector<Vector3>   A;
        std::vector<float>   b;
        foreach3D(Vector3i(vSize.x, vSize.y, vSize.z), [&](const Vector3i& p) {
            std::vector<bool> inside(vCorners.size(), false);
            int  numInside = 0;
            for (unsigned int ci = 0; ci < vCorners.size(); ++ci) {
                inside[ci] = (buffer.get_voxel_f(p + vCorners[ci], VoxelBuffer::CHANNEL_SDF) < 0);
                if (inside[ci]) {
                    numInside += 1;
                }
            }

            if (numInside == 0 || numInside == vCorners.size()) {
                return;
            }

            std::vector<bool> crossingCorners(vCorners.size(), false);

            for (int ai = 0; ai < NumEdges; ++ai) {
                auto&& e = vEdges[ai];
                if (inside[e[0]] != inside[e[1]]) {
                    crossingCorners[e[0]] = true;
                    crossingCorners[e[1]] = true;
                }
            }

            planes.clear();
            A.clear();
            b.clear();

            // Add candidate planes from all corners with a sign-changing edge:
            std::vector<Vector3> vNormals = getNormal(buffer, p);
            for (unsigned int ci = 0; ci < vCorners.size(); ++ci) {
                if (!crossingCorners[ci]) {
                    // Corner has no edge with sign change
                    continue;
                }

                Vector3i p_ni = p + vCorners[ci];
                const Vector3 p_n = Vector3(p_ni.x, p_ni.y, p_ni.z);
                float fValue = buffer.get_voxel_f(p_n, VoxelBuffer::CHANNEL_SDF);
                Vector3 vNormal = vNormals[ci];
                vNormal.normalize();
                if (vNormal == Vector3(0, 0, 0)) {
                    continue;
                }

                Plane plane = { fValue,  vNormal };

                if (std::fabs(plane.dist) > MaxCornerDist) {
                    // Large distance change - produced by bad input
                    // This will most likely cause a bad vertex, so skip this corner
                    continue;
                }

                plane.dist = plane.normal.dot(p_n) - plane.dist;
                planes.push_back(plane);
            }

            for (int ai = 0; ai < 3; ++ai) {
                Vector3 normal = (CenterPush * vAxes[ai]).to_vec3();
                Vector3 pos = p.to_vec3() + Vector3(0.5, 0.5, 0.5);
                planes.push_back(Plane{ normal.dot(pos), normal });
            }

            for (auto&& p : planes) {
                A.push_back(p.normal);
                b.push_back(p.dist);
            }

            Vector3 vertex = LeastSquareSolver::leastSquares(A.size(), A.data(), b.data());

            auto voxelCenter = p.to_vec3() + Vector3(0.5f, 0.5f, 0.5f);
            // leastSquares failed
            if (!LeastSquareSolver::isFinite(vertex)) {
                float count = 1.0f;
                Vector3 vVertexModify = voxelCenter;
                for (unsigned int x = 0; x < vCorners.size(); x++) {
                    Vector3 p_n = (p + vCorners[x]).to_vec3();
                    if (buffer.get_voxel_f(p_n, VoxelBuffer::CHANNEL_SDF) <= 0) {
                        vVertexModify += p_n;
                        count += 1.0f;
                    }
                }
                vertex = vVertexModify / count;
            }

            Vector3 clamped = clamp_to(vertex, p, p + Vector3i(1, 1, 1));
            if (bClamp)
            {
                if (vertex != clamped) {
                    vertex = Vector3i(clamped.x, clamped.y, clamped.z).to_vec3();
                }
            }
            else if (voxelCenter.distanceTo(vertex) > FarAway) {
                float count = 1.0f;
                Vector3 vVertexModify = voxelCenter;
                for (unsigned int x = 0; x < vCorners.size(); x++) {
                    Vector3 p_n = (p + vCorners[x]).to_vec3();
                    if (buffer.get_voxel_f(p_n, VoxelBuffer::CHANNEL_SDF) <= 0) {
                        vVertexModify += p_n;
                        count += 1.0f;
                    }
                }
                vertex = vVertexModify / count;
            }
            //minutes min_padding
            vertex -= Vector3(min_padding, min_padding, min_padding);

            //vertIdx.set_voxel(mesh.vecs.size(), p.x, p.y, p.z, VoxelBuffer::CHANNEL_DATA3);
            mesh.vecs.push_back(vertex);
            mesh.normals.push_back(Vector3(0));
            if (p.x == 0 || p.y == 0 || p.z == 0) {
                mesh.valid.push_back(false);
            }
            else if (p.x == blocksize_with_padding - 2 || p.y == blocksize_with_padding - 2 || p.z == blocksize_with_padding - 2) {
                mesh.valid.push_back(false);
            }
            else {
                mesh.valid.push_back(true);
            }
        });
    }

    void DualContouring::constructFaces(const VoxelBuffer &buffer, TriMesh &mesh, VoxelBuffer &vertIdx, const int blocksize_with_padding) {
        Vector3i vSize = buffer.get_size();
        std::vector<std::vector<Vector3i>> vNeighbors = {
                    {{0, 0, 1}, {0, 1, 0}, {0, 1, 1}},
                    {{0, 0, 1}, {1, 0, 0}, {1, 0, 1}},
                    {{0, 1, 0}, {1, 0, 0}, {1, 1, 0}}
        };
        foreach3D(Vector3i(vSize.x, vSize.y, vSize.z), [&](const Vector3i &p) {
            uint64_t v0 = vertIdx.get_voxel(Vector3i(p.x, p.y, p.z), VoxelBuffer::CHANNEL_DATA3);
            if (v0 == MAX_UNSIGNED_LONG) {
                return;
            }

            std::vector<bool> inside((int)vCorners.size(), false);

            for (unsigned int ci = 0; ci < vCorners.size(); ++ci) {
                inside[ci] = (buffer.get_voxel_f(p + vCorners[ci], VoxelBuffer::CHANNEL_SDF) < 0);
            }

            for (int ai = 0; ai < 3; ++ai) {
                auto&& e = FarEdges[ai];
                if (inside[e[0]] == inside[e[1]]) {
                    continue;  // Not a crossing
                }

                uint64_t v[3];

                for (int j = 0; j < 3; j++) {
                    Vector3i np = p + vNeighbors[ai][j];
                    Vector3i vp = Vector3i(np.x, np.y, np.z);

                    v[j] = vertIdx.get_voxel(vp, VoxelBuffer::CHANNEL_DATA3);
                }

                if (v[0] == MAX_UNSIGNED_LONG || v[1] == MAX_UNSIGNED_LONG || v[2] == MAX_UNSIGNED_LONG) {
                    continue;
                }

                Vector3i t0 = { (int)v0, (int)v[0], (int)v[2] };
                Vector3i t1 = { (int)v0, (int)v[2], (int)v[1] };

                // Get the normals right:
                if (inside[e[0]] != (ai == 1)) {
                    t0 = flip(t0);
                    t1 = flip(t1);
                }

                Vector3 vNormal = Triangle(mesh.vecs[t0.x], mesh.vecs[t0.y], mesh.vecs[t0.z]).getNormal();
                mesh.normals[t0.x] += vNormal;
                mesh.normals[t0.x].normalize();
                mesh.normals[t0.y] += vNormal;
                mesh.normals[t0.y].normalize();
                mesh.normals[t0.z] += vNormal;
                mesh.normals[t0.z].normalize();
                vNormal = Triangle(mesh.vecs[t1.x], mesh.vecs[t1.y], mesh.vecs[t1.z]).getNormal();
                mesh.normals[t1.x] += vNormal;
                mesh.normals[t1.x].normalize();
                mesh.normals[t1.y] += vNormal;
                mesh.normals[t1.y].normalize();
                mesh.normals[t1.z] += vNormal;
                mesh.normals[t1.z].normalize();

                if (mesh.valid[t0.x] && mesh.valid[t0.y] && mesh.valid[t0.z]) {
                    mesh.triangles.push_back(flip(t0));
                }
                if (mesh.valid[t1.x] && mesh.valid[t1.y] && mesh.valid[t1.z]) {
                    mesh.triangles.push_back(flip(t1));
                }
            }
        });
    }

    TriMesh DualContouring::dualContouring(const VoxelBuffer &buffer,
                                           const int min_padding, const int max_padding,
                                           const int blocksize_with_padding) {
        TriMesh mesh;
//        VoxelBuffer vertIdx;
//        vertIdx.create(buffer.get_size());

//        vertIdx.set_channel_depth(VoxelBuffer::CHANNEL_DATA3, VoxelBuffer::DEPTH_64_BIT);
//        vertIdx.fill(MAX_UNSIGNED_LONG, VoxelBuffer::CHANNEL_DATA3);

//        constructVertices(buffer, mesh, vertIdx, min_padding, blocksize_with_padding);
//        constructFaces(buffer, mesh, vertIdx, blocksize_with_padding);

        return mesh;
    }
};

