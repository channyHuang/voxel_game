#ifndef SURFACENETS_H
#define SURFACENETS_H

#include <vector>
#include <functional>
#include <unordered_map>
#include <memory>
#include <fstream>

#include <QDebug>

#ifndef Real
#define Real float
#endif

#include "../common_math/vector3.h"
#include "../common_math/vector3i.h"
#include "../voxels/voxel_buffer.h"
#include "common_math/triangle.h"

using namespace std;

namespace SURFACE_NETS {

    enum MaterialType {
        AIR = 0,
        ROAD,
        STONE,
        GRASSLAND,
        GRASS,
        SAND,
        LAVA,
        WATER,
        BRICK,
        WOOD
    };

	enum MeshType {
		MESH_UNKNOWN = -1,
		MESH_SMOOTH = 0,
		MESH_HALFSMOOTH,
		MESH_CUBE,
		MESH_WATER,
		MESH_QUANTIZE
	};
	struct VertexMesh {
		vector<Vector3> vertices_;
		vector<Vector3i> faces_; //vertices index
		vector<Vector3> normals_;
		vector<uint8_t> vertexValid_; //vertices in border are not valid
		vector<Vector3> face_normals_;
		vector<uint8_t> share_point_;
		vector<Vector3i> edge_faces_;
		//mesh with water
		vector<MeshType> mesh_type_;
		vector<Vector3i> solid_faces_; //non-water faces
		vector<Vector3i> water_faces_; //water faces

		void clear() {
			vertices_.clear();
			faces_.clear();
			normals_.clear();
			vertexValid_.clear();
			face_normals_.clear();
			share_point_.clear();
			edge_faces_.clear();

			mesh_type_.clear();
			solid_faces_.clear();
			water_faces_.clear();
		}

		void swap(VertexMesh& mesh) {
			this->clear();
			vertices_.swap(mesh.vertices_);
			faces_.swap(mesh.faces_);
			normals_.swap(mesh.normals_);
			vertexValid_.swap(mesh.vertexValid_);
			face_normals_.swap(mesh.face_normals_);
			share_point_.swap(mesh.share_point_);
			edge_faces_.swap(mesh.edge_faces_);

			mesh_type_.swap(mesh.mesh_type_);
			solid_faces_.swap(mesh.solid_faces_);
			water_faces_.swap(mesh.water_faces_);
		}

		void reserve(const size_t capacity) {
			vertices_.reserve(capacity);
			faces_.reserve(capacity);
			normals_.reserve(capacity);
			vertexValid_.reserve(capacity);
			face_normals_.reserve(capacity);
			share_point_.reserve(capacity);
			edge_faces_.reserve(capacity);

			mesh_type_.reserve(capacity);
			solid_faces_.reserve(capacity);
			water_faces_.reserve(capacity);
		}

		bool isWater(Vector3i face_index) {
			return (mesh_type_[face_index.x] == MeshType::MESH_WATER
				|| mesh_type_[face_index.y] == MeshType::MESH_WATER
				|| mesh_type_[face_index.z] == MeshType::MESH_WATER);
		}

        void outputToObjFile(std::string sFileName, Vector3i position, bool bseperate = false) {
            if (vertices_.size() <= 0) return;
            Vector3i offset = position * 16 - Vector3i(2);
            if (!bseperate) {
                if (faces_.size() <= 0) return;

                std::ofstream ofs(sFileName + std::to_string(position.x) + "_" + std::to_string(position.y) + "_" + std::to_string(position.z) + ".obj");
                for (size_t i = 0; i < vertices_.size(); i++) {
                    ofs << "v " << vertices_[i][0] + offset.x << " " << vertices_[i][1] + offset.y << " " << vertices_[i][2] + offset.z << std::endl;
                    ofs << "vn " << normals_[i].x << " " << normals_[i].y << " " << normals_[i].z << std::endl;
                }
                for (size_t i = 0; i < faces_.size(); ++i) {
                    ofs << "f " << faces_[i].x + 1 << " " << faces_[i].y + 1 << " " << faces_[i].z + 1 << std::endl;
                }
                ofs.close();
            } else {
                if (water_faces_.size() <= 0 && solid_faces_.size() <= 0) return;
                if (water_faces_.size() > 0) {
                    std::ofstream ofs_water(sFileName + "_water_" + std::to_string(position.x) + "_" + std::to_string(position.y) + "_" + std::to_string(position.z) + ".obj");
                    for (size_t i = 0; i < vertices_.size(); i++) {
                        ofs_water << "v " << vertices_[i][0] + offset.x << " " << vertices_[i][1] + offset.y << " " << vertices_[i][2] + offset.z << std::endl;
                        ofs_water << "vn " << normals_[i].x << " " << normals_[i].y << " " << normals_[i].z << std::endl;
                    }
                    for (size_t i = 0; i < water_faces_.size(); ++i) {
                        ofs_water << "f " << water_faces_[i].x + 1 << " " << water_faces_[i].y + 1 << " " << water_faces_[i].z + 1 << std::endl;
                    }
                    ofs_water.close();
                }
                if (solid_faces_.size() > 0) {
                    std::ofstream ofs_water(sFileName + "_solid_" + std::to_string(offset.x) + "_" + std::to_string(offset.y) + "_" + std::to_string(offset.z) + ".obj");
                    for (size_t i = 0; i < vertices_.size(); i++) {
                        ofs_water << "v " << vertices_[i][0] + offset.x << " " << vertices_[i][1] + offset.y << " " << vertices_[i][2] + offset.z << std::endl;
                        ofs_water << "vn " << normals_[i].x << " " << normals_[i].y << " " << normals_[i].z << std::endl;
                    }
                    for (size_t i = 0; i < solid_faces_.size(); ++i) {
                        ofs_water << "f " << solid_faces_[i].x + 1 << " " << solid_faces_[i].y + 1 << " " << solid_faces_[i].z + 1 << std::endl;
                    }
                    ofs_water.close();
                }
            }
        }
	};

	namespace {
		static const unsigned char kVertexIndexTable[][3] = {
			{0, 0, 0},
			{1, 0, 0},
			{1, 1, 0},
			{0, 1, 0},

			{0, 0, 1},
			{1, 0, 1},
			{1, 1, 1},
			{0, 1, 1},
		};

		static const unsigned char kEdgeVertexTable[12][2][3] = {
			{ {0, 0, 0}, {1, 0, 0} },
			{ {1, 0, 0}, {1, 1, 0} },
			{ {1, 1, 0}, {0, 1, 0} },
			{ {0, 1, 0}, {0, 0, 0} },

			{ {0, 0, 1}, {1, 0, 1} },
			{ {1, 0, 1}, {1, 1, 1} },
			{ {1, 1, 1}, {0, 1, 1} },
			{ {0, 1, 1}, {0, 0, 1} },

			{ {0, 0, 0}, {0, 0, 1} },
			{ {1, 0, 0}, {1, 0, 1} },
			{ {1, 1, 0}, {1, 1, 1} },
			{ {0, 1, 0}, {0, 1, 1} },
		};

		// Note: indices here are arranged in the order "ABC CBD" ("strip order") to make it easy to extract diagonal later
		static const unsigned int kOffsetTable[2][2][6] = {
			{
				{ 1, 0, 2, 2, 0, 3 },
				{ 1, 2, 0, 0, 2, 3 },
			},
			{
				{ 0, 3, 1, 1, 3, 2 },
				{ 0, 1, 3, 3, 1, 2 },
			}
		};

		//6561 = 3^8
		static unsigned short gEdgeTable[6561];
	}

	bool gen_vertex_by_voxel(
		std::function<Real(const Vector3&)> const& sdfFunction,
        std::function<MaterialType(const Vector3&)> const& materialFunction,
        const Vector3i &vpos,
        Vector3 &mesh_vertex,
		bool &bSharePoint,
		const int lod = 0,
		float const isovalue = 0.f);

	bool gen_avg_vertex_by_voxel(
		std::function<Real(const Vector3&)> const& sdfFunction,
        std::function<MaterialType(const Vector3&)> const& materialFunction,
        const Vector3i &vpos,
        Vector3 &mesh_vertex,
		bool &bSharePoint,
		const int lod = 0,
		float const isovalue = 0.f);

    std::shared_ptr<VertexMesh> surface_nets(
        std::function<Real(const Vector3i&)> const& sdfFunction,
        std::function<MaterialType(const Vector3&)> const& materialFunction,
        const Vector3i &vVoxelSize,
		const int nMinPadding, const int nMaxPadding,
		const int lod = 0,
        const Vector3i &start_position = Vector3i(0),
		float const isovalue = 0.f);

	std::shared_ptr<VertexMesh> surface_nets_reduce_surface(
        std::function<Real(const Vector3i&)> const& sdfFunction,
        std::function<MaterialType(const Vector3&)> const& materialFunction,
        const Vector3i &vVoxelSize,
		const int nMinPadding, const int nMaxPadding,
		const int lod = 0,
        const Vector3i &position = Vector3i(0),
		float const isovalue = 0.f);

	std::shared_ptr<VertexMesh> surface_nets_lod_from_upper_lod(
		std::shared_ptr<VertexMesh> pmesh,
		std::function<Real(const Vector3&)> const& sdfFunction,
        std::function<MaterialType(const Vector3&)> const& materialFunction,
        const Vector3i &vVoxelSize,
		const int nMinPadding, const int nMaxPadding,
		int downscale_lod);
	//position of vertex in voxels when material type is regular material
    static Real fcube_position = 0.8f;


	//with water
	void prepareTables();

    inline static Vector3 vector3i2Vector3(const Vector3i &vposi) {
        return Vector3(static_cast<Real>(vposi.x), static_cast<Real>(vposi.y), static_cast<Real>(vposi.z));
	}
	inline static Vector3 round(const Vector3& v) {
		Real x = static_cast<Real>(v.x < 0 ? std::ceil(v.x - 0.5) : std::floor(v.x + 0.5));
		Real y = static_cast<Real>(v.y < 0 ? std::ceil(v.y - 0.5) : std::floor(v.y + 0.5));
		Real z = static_cast<Real>(v.z < 0 ? std::ceil(v.z - 0.5) : std::floor(v.z + 0.5));
		return Vector3(x, y, z);
	}
	
	inline void hash_combine(std::size_t& seed) { }
	template <typename T, typename... Rest>
	inline void hash_combine(std::size_t& seed, const T& vposi, Rest... rest) {
		size_t hash_v = vposi.y & 255 | (vposi.x & 32767) << 8 | (vposi.z & 32767) << 24;
		seed ^= hash_v + 0x9e3779b9 + (seed << 6) + (seed >> 2);
		hash_combine(seed, rest...);
	}

	static size_t avalanche(size_t v) {
		v += ~(v << 15);
		v ^= (v >> 10);
		v += (v << 3);
		v ^= (v >> 6);
		v += ~(v << 11);
		v ^= (v >> 16);
		return v;
	}

	static size_t computeSeed(const Vector3i& vposi) {
		size_t result = 0;
		hash_combine(result, vposi);
		return avalanche(result);
	}

    static Vector3i vector3FloorOrCeil(const Vector3 &pos, bool bfloor = true) {
		if (bfloor) {
            Vector3 &&pos2floor = pos.getFloor();
            Vector3i &&vposi = Vector3i(static_cast<int>(pos2floor.x), static_cast<int>(pos2floor.y), static_cast<int>(pos2floor.z));
			return vposi;
		}
		//bceil
        Vector3 &&pos2ceil = Vector3(pos);
		pos2ceil.ceil();
        Vector3i &&vposi = Vector3i(static_cast<int>(pos2ceil.x), static_cast<int>(pos2ceil.y), static_cast<int>(pos2ceil.z));
		return vposi;
	}

	static Vector3 computePoint(const Vector3 &vsmooth, float cellSize, MeshType material, size_t seed);

    inline std::pair<size_t, MaterialType> reduceMaterials(
        const std::pair<size_t, MaterialType>& m0,
        const std::pair<size_t, MaterialType>& m1) {
		if (m0.second == m1.second)
			return std::make_pair(m0.first + m1.first, m0.second);
		else if (m0.first != m1.first)
			return m0.first > m1.first ? m0 : m1;
		else
			return m0.second < m1.second ? m0 : m1;
	}

	static void pushQuad(VertexMesh &mesh,
		size_t v0, size_t v1, size_t v2, size_t v3,
		bool flip) {
		size_t v[] = { v1, v2, v3, v0 };
		MeshType m[4];
		for (int i = 0; i < 4; ++i) {
			m[i] = mesh.mesh_type_[v[i]];
		}
		bool bflip = (m[1] == m[3]) && (m[1] == m[2] || m[1] == m[0]);
		const unsigned int* offset = kOffsetTable[bflip][flip];

		Vector3i vface_index = Vector3i(v[offset[0]], v[offset[1]], v[offset[2]]);
		mesh.faces_.push_back(vface_index);
		vface_index = Vector3i(v[offset[3]], v[offset[4]], v[offset[5]]);
		mesh.faces_.push_back(vface_index);
	}

    void separateSolidAndWaterFaces(std::shared_ptr<VertexMesh> pmesh);
    std::shared_ptr<VertexMesh> surface_nets_with_water(
		std::function<Real(const Vector3&)> const& sdfFunction,
        std::function<MaterialType(const Vector3&)> const& materialFunction,
        std::function<std::pair<int, int>(const Vector3i&)> const& funVoxelOccupancyAndTag,
		std::function<int(const Vector3i&)> const& funVoxelTag,
        const Vector3i &vVoxelSize,
		const int lod);

	static std::pair<Vector3, MeshType> gen_avg_vertex_by_voxel(
		std::function<Real(const Vector3&)> const& sdfFunction,
        std::function<MaterialType(const Vector3&)> const& materialFunction,
        std::function<std::pair<int, int>(const Vector3i&)> const& funVoxelOccupancyAndTag,
        const Vector3i &vVoxelSize,
		const int lod,
		const Vector3i &vposi,
		const int edgemask);
	static void gen_faces(VertexMesh &mesh, 
        const Vector3i &vvoxel_size,
        std::function<std::pair<int, int>(const Vector3i&)> const& funVoxelOccupancyAndTag,
		std::unordered_map<std::size_t, std::size_t> &mvoxelIndex2MeshIndex);
}

#endif
