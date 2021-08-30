#include "surface_nets.h"

#include <QDebug>

/*
*
* Visit every voxel in the regular vsize with the following configuration
*
* Coordinate frame
*
*       z y
*       |/
*       o--x
*
* Voxel corner indices
*
*        7          6
*        o----------o
*       /|         /|
*     4/ |       5/ |
*     o--|-------o  |
*     |  o-------|--o
*     | /3       | /2
*     |/         |/
*     o----------o
*     0          1
*
*
* we define a  lower boundary cube to be a cube with (i=0       || j=0       || k=0	    )
* we define an upper boundary cube to be a cube with (i >= sx-1 || j >= sy-1 || k >=sz-1)
*
* lower boundary cubes have missing neighbor voxels, so we don't triangulate
* when the current voxel is a boundary cube. Our method of quad generation considers
* the following potentially crossed edges of an active cube:
*
* Potentially crossed edges
*
*     4
*     o
*     |  o
*     | /3
*     |/
*     o----------o
*     0          1
*
* For potentially crossed edges, we generate quads connecting the 4 cubes of
* the crossed edge, so the considered neighbors are:
*
*
*           o----------o----------o
*          /|         /|         /|
*         / |   0    / |        / | <---- this is the current active cube (i,j,k) in the above coordinate frame
*        o--|-------o--|-------o  |
*       /|  o------/|--o------/|--o
*      / | /|1    / | /|2    / | /|
*     o--|/-|--5-o--|/-|--4-o  |/ |
*     |  o--|----|--o--|----|--o  |
*     | /|  o----|-/|--o----|-/|--o
*     |/ | /     |/ | /3    |/ | /
*     o--|/------o--|/------o  |/
*        o-------|--o-------|--o
*                | /        | /
*                |/         |/
*                o----------o
*
*/



TriMesh SurfaceNets::surfaceNets(const std::function<float (float, float, float)> &sdfFunction,
                     Vector3i vsize,
                     float const isovalue)
{
    TriMesh mesh{};

    // mapping from 3d coordinates to 1d
    auto const get_active_cube_index =
        [](std::size_t x, std::size_t y, std::size_t z, Vector3i const& vsize) -> std::size_t {
        return x + (y * vsize.x) + (z * vsize.x * vsize.y);
    };

    auto const get_ijk_from_idx =
        [](std::size_t active_cube_index,
           Vector3i const& vsize) -> std::tuple<std::size_t, std::size_t, std::size_t> {
        std::size_t i = (active_cube_index) % vsize.x;
        std::size_t j = (active_cube_index / vsize.x) % vsize.y;
        std::size_t k = (active_cube_index) / (vsize.x * vsize.y);
        return std::make_tuple(i, j, k);
    };

    auto const is_scalar_positive = [](float scalar, float isovalue) -> bool
    {
        return scalar >= isovalue;
    };

    auto const are_edge_scalars_bipolar = [&is_scalar_positive](float scalar1, float scalar2, float isovalue) -> bool
    {
        return is_scalar_positive(scalar1, isovalue) != is_scalar_positive(scalar2, isovalue);
    };

    auto const is_lower_boundary_cube = [](std::size_t i, std::size_t j, std::size_t k)
    {
        return (i == 0 || j == 0 || k == 0);
    };


    // mapping from active cube indices to vertex indices of the generated mesh
    std::unordered_map<std::size_t, std::uint64_t> active_cube_to_vertex_index_map{};

    for (int k = 0; k < vsize.z - 1; ++k) {
        for (int j = 0; j < vsize.y - 1; ++j) {
            for (int i = 0; i < vsize.x - 1; ++i) {
                Vector3 const voxel_corner_vsize_positions[8] =
                {
                    { static_cast<float>(i)    , static_cast<float>(j)    , static_cast<float>(k)     },
                    { static_cast<float>(i + 1), static_cast<float>(j)    , static_cast<float>(k)     },
                    { static_cast<float>(i + 1), static_cast<float>(j + 1), static_cast<float>(k)     },
                    { static_cast<float>(i)    , static_cast<float>(j + 1), static_cast<float>(k)     },
                    { static_cast<float>(i)    , static_cast<float>(j)    , static_cast<float>(k + 1) },
                    { static_cast<float>(i + 1), static_cast<float>(j)    , static_cast<float>(k + 1) },
                    { static_cast<float>(i + 1), static_cast<float>(j + 1), static_cast<float>(k + 1) },
                    { static_cast<float>(i)    , static_cast<float>(j + 1), static_cast<float>(k + 1) },
                };

                Vector3 voxel_corner_positions[8];
                for (int x = 0; x < 8; x++) {
                    voxel_corner_positions[x] = voxel_corner_vsize_positions[x];
                }

                float const voxel_corner_values[8] =
                {
                    sdfFunction(voxel_corner_positions[0].x, voxel_corner_positions[0].y, voxel_corner_positions[0].z),
                    sdfFunction(voxel_corner_positions[1].x, voxel_corner_positions[1].y, voxel_corner_positions[1].z),
                    sdfFunction(voxel_corner_positions[2].x, voxel_corner_positions[2].y, voxel_corner_positions[2].z),
                    sdfFunction(voxel_corner_positions[3].x, voxel_corner_positions[3].y, voxel_corner_positions[3].z),
                    sdfFunction(voxel_corner_positions[4].x, voxel_corner_positions[4].y, voxel_corner_positions[4].z),
                    sdfFunction(voxel_corner_positions[5].x, voxel_corner_positions[5].y, voxel_corner_positions[5].z),
                    sdfFunction(voxel_corner_positions[6].x, voxel_corner_positions[6].y, voxel_corner_positions[6].z),
                    sdfFunction(voxel_corner_positions[7].x, voxel_corner_positions[7].y, voxel_corner_positions[7].z)
                };

                bool const edge_bipolarity_array[12] =
                {
                    are_edge_scalars_bipolar(voxel_corner_values[edges[0][0]], voxel_corner_values[edges[0][1]], isovalue),
                    are_edge_scalars_bipolar(voxel_corner_values[edges[1][0]], voxel_corner_values[edges[1][1]], isovalue),
                    are_edge_scalars_bipolar(voxel_corner_values[edges[2][0]], voxel_corner_values[edges[2][1]], isovalue),
                    are_edge_scalars_bipolar(voxel_corner_values[edges[3][0]], voxel_corner_values[edges[3][1]], isovalue),
                    are_edge_scalars_bipolar(voxel_corner_values[edges[4][0]], voxel_corner_values[edges[4][1]], isovalue),
                    are_edge_scalars_bipolar(voxel_corner_values[edges[5][0]], voxel_corner_values[edges[5][1]], isovalue),
                    are_edge_scalars_bipolar(voxel_corner_values[edges[6][0]], voxel_corner_values[edges[6][1]], isovalue),
                    are_edge_scalars_bipolar(voxel_corner_values[edges[7][0]], voxel_corner_values[edges[7][1]], isovalue),
                    are_edge_scalars_bipolar(voxel_corner_values[edges[8][0]], voxel_corner_values[edges[8][1]], isovalue),
                    are_edge_scalars_bipolar(voxel_corner_values[edges[9][0]], voxel_corner_values[edges[9][1]], isovalue),
                    are_edge_scalars_bipolar(voxel_corner_values[edges[10][0]], voxel_corner_values[edges[10][1]], isovalue),
                    are_edge_scalars_bipolar(voxel_corner_values[edges[11][0]], voxel_corner_values[edges[11][1]], isovalue),
                };

                // an active voxel must have at least one bipolar edge
                bool const is_voxel_active = edge_bipolarity_array[0] ||
                    edge_bipolarity_array[1] ||
                    edge_bipolarity_array[2] ||
                    edge_bipolarity_array[3] ||
                    edge_bipolarity_array[4] ||
                    edge_bipolarity_array[5] ||
                    edge_bipolarity_array[6] ||
                    edge_bipolarity_array[7] ||
                    edge_bipolarity_array[8] ||
                    edge_bipolarity_array[9] ||
                    edge_bipolarity_array[10] ||
                    edge_bipolarity_array[11];

                // cubes that are not active do not generate mesh vertices
                if (!is_voxel_active)
                    continue;

                // store all edge intersection points with the implicit surface in voxel vsize coordinates
                std::vector<Vector3> edge_intersection_points;

                // visit every bipolar edge
                for (std::size_t e = 0; e < 12; ++e)
                {
                    if (!edge_bipolarity_array[e])
                        continue;

                    // get points p1, p2 of the edge e in vsize coordinates
                    auto p1 = voxel_corner_vsize_positions[edges[e][0]];
                    auto p2 = voxel_corner_vsize_positions[edges[e][1]];

                    // get value of the implicit function at edge vertices
                    auto s1 = voxel_corner_values[edges[e][0]];
                    auto s2 = voxel_corner_values[edges[e][1]];

                    // perform linear interpolation using implicit function
                    // values at vertices
                    float t = (isovalue - s1) / (s2 - s1);
                    edge_intersection_points.emplace_back(p1 + t * (p2 - p1));
                }

                float number_of_intersection_points = static_cast<float>(edge_intersection_points.size());
                Vector3 sum_of_intersection_points = std::accumulate(
                    edge_intersection_points.cbegin(),
                    edge_intersection_points.cend(),
                    Vector3{ 0.f, 0.f, 0.f }
                );
                Vector3 geometric_center_of_edge_intersection_points = sum_of_intersection_points / number_of_intersection_points;
                Vector3 const mesh_vertex = geometric_center_of_edge_intersection_points;

                std::size_t const active_cube_index = get_active_cube_index(i, j, k, vsize);
                std::uint64_t const vertex_index = mesh.vertices.size();
                active_cube_to_vertex_index_map[active_cube_index] = vertex_index;

                mesh.vertices.emplace_back(mesh_vertex);
            }
        }
    }

    // Triangulation
	for (auto const& key_value : active_cube_to_vertex_index_map)
	{
		std::size_t   const active_cube_index = key_value.first;
		std::uint64_t const vertex_index = key_value.second;

        auto const ijk = get_ijk_from_idx(active_cube_index, vsize);
        int i = std::get<0>(ijk);
        int j = std::get<1>(ijk);
        int k = std::get<2>(ijk);

        if (is_lower_boundary_cube(i, j, k))
			continue;

        int const neighbor_vsize_positions[6][3] =
		{
			{ i - 1, j    , k     },
			{ i - 1, j - 1, k     },
			{ i    , j - 1, k     },
			{ i    , j - 1, k - 1 },
			{ i    , j    , k - 1 },
			{ i - 1, j    , k - 1 }
		};

        Vector3i const voxel_corners_of_interest[4] =
		{
            // vertex 0, 4, 3, 1
            {i, j, k},
            {i, j, k + 1},
            {i, j + 1, k},
            {i + 1, j, k}
        };

		float const edge_scalar_values[3][2] =
		{
			// directed edge (0,4)
			{
                sdfFunction(voxel_corners_of_interest[0].x, voxel_corners_of_interest[0].y, voxel_corners_of_interest[0].z),
                sdfFunction(voxel_corners_of_interest[1].x, voxel_corners_of_interest[1].y, voxel_corners_of_interest[1].z)
			},
			// directed edge (3,0)
			{
                sdfFunction(voxel_corners_of_interest[2].x, voxel_corners_of_interest[2].y, voxel_corners_of_interest[2].z),
                sdfFunction(voxel_corners_of_interest[0].x, voxel_corners_of_interest[0].y, voxel_corners_of_interest[0].z)
			},
			// directed edge (0,1)
			{
                sdfFunction(voxel_corners_of_interest[0].x, voxel_corners_of_interest[0].y, voxel_corners_of_interest[0].z),
                sdfFunction(voxel_corners_of_interest[3].x, voxel_corners_of_interest[3].y, voxel_corners_of_interest[3].z)
			}
		};



		// look at each potentially generated quad
		for (std::size_t i = 0; i < 3; ++i)
		{
			auto const neighbor1 = get_active_cube_index(
                neighbor_vsize_positions[quad_neighbors[i][0]][0],
                neighbor_vsize_positions[quad_neighbors[i][0]][1],
                neighbor_vsize_positions[quad_neighbors[i][0]][2],
                vsize);

			auto const neighbor2 = get_active_cube_index(
                neighbor_vsize_positions[quad_neighbors[i][1]][0],
                neighbor_vsize_positions[quad_neighbors[i][1]][1],
                neighbor_vsize_positions[quad_neighbors[i][1]][2],
                vsize);

			auto const neighbor3 = get_active_cube_index(
                neighbor_vsize_positions[quad_neighbors[i][2]][0],
                neighbor_vsize_positions[quad_neighbors[i][2]][1],
                neighbor_vsize_positions[quad_neighbors[i][2]][2],
                vsize);

			if (active_cube_to_vertex_index_map.count(neighbor1) == 0 ||
				active_cube_to_vertex_index_map.count(neighbor2) == 0 ||
				active_cube_to_vertex_index_map.count(neighbor3) == 0)
				continue;

			std::size_t const neighbor_vertices[3] =
			{
				active_cube_to_vertex_index_map.at(neighbor1),
				active_cube_to_vertex_index_map.at(neighbor2),
				active_cube_to_vertex_index_map.at(neighbor3)
			};

			auto const& neighbor_vertices_order =
				edge_scalar_values[i][1] > edge_scalar_values[i][0] ?
				quad_neighbor_orders[0] :
				quad_neighbor_orders[1];

            int v0 = vertex_index;
            int v1 = neighbor_vertices[neighbor_vertices_order[0]];
            int v2 = neighbor_vertices[neighbor_vertices_order[1]];
            int v3 = neighbor_vertices[neighbor_vertices_order[2]];

			mesh.add_face({ v0, v1, v2 });
			mesh.add_face({ v0, v2, v3 });
		}
	}

    return mesh;
}

