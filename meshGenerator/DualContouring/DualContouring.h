#pragma once

#include "vector2i.h"
#include "vector3i.h"
#include "triangle.h"
#include <iostream>
#include <vector>
#include <array>
#include "LeastSquareSolver.h"
#include "voxel_buffer.h"
//using MC Tables to calculate normals quickly
#include "TransvoxelTables.hpp"

namespace DualContouring {
	struct TriMesh {
		std::vector<Vector3>  vecs;
		std::vector<Vector3i>  triangles;
		std::vector<Vector3> normals;
		std::vector<bool> valid;
	};

	class DualContouring {
	public:
		const bool bClamp = false;
		const real MaxCornerDist = 1.0;
		const real FarAway = 1.0f; // 2.5;
		const real CenterPush = 0.01;

		const int NumCorners = 8;
		const int NumEdges = 12;
		const unsigned int MAX_UNSIGNED_LONG = 4294967295U;

		const std::vector<Vector3i> vCorners = {
			{0,0,0}, {0,0,1}, {0,1,0}, {0,1,1},
			{1,0,0}, {1,0,1}, {1,1,0}, {1,1,1},
		};

        const Vector3i vAxes[3] = {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}};
		// Indices into Corners:
        std::vector<Vector2i> vEdges = {
			{0,1}, {0,2}, {0,4},
			{1,3}, {1,5},
			{2,3}, {2,6},
			{3,7},
			{4,5}, {4,6},
			{5,7},
			{6,7}
		};

		struct Plane {
			real dist;
			Vector3 normal;
		};

		// The edges leading to our far corner:
        const Vector2i FarEdges[3] = { {3,7}, {5,7}, {6,7} };

		inline Vector3i flip(Vector3i p) {
            return Vector3i(p[0], p[2], p[1]);
		}

		std::vector<Vector3> getNormal(const VoxelBuffer &voxels, Vector3i pos);
		TriMesh dualContouring(const VoxelBuffer &buffer, const int min_padding, const int max_padding, const int blocksize_with_padding);
		void constructFaces(const VoxelBuffer &buffer, TriMesh &mesh, VoxelBuffer &vertIdx, const int blocksize_with_padding);
		void constructVertices(const VoxelBuffer &buffer, TriMesh &mesh, VoxelBuffer &vertIdx, const int min_padding, const int blocksize_with_padding);

		template<typename Fun>
        void foreach2D(Vector2i min,
            Vector2i max,
			const Fun& fun)
		{
			for (int y = min.y; y < max.y; ++y)
				for (int x = min.x; x < max.x; ++x)
                    fun(Vector2i(x, y));
		}

		template<typename Fun>
        void foreach2D(Vector2i max, const Fun& fun) {
			foreach2D({ 0,0 }, max, fun);
		}

		template<typename Fun>
		void foreach3D(Vector3i min,
			Vector3i max,
			const Fun& fun)
		{
			for (int z = min.z; z < max.z; ++z)
				for (int y = min.y; y < max.y; ++y)
					for (int x = min.x; x < max.x; ++x)
						fun(Vector3i(x, y, z));
		}

		template<typename Fun>
		void foreach3D(Vector3i max, const Fun& fun) {
			foreach3D({ 0,0,0 }, max, fun);
		}
	};
}
