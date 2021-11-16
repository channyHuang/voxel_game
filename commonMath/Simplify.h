#pragma once

#include "vector3.h"
#include <vector>

namespace Simplify {
	class SymetricMatrix {
	public:
        SymetricMatrix(float c = 0) { for (unsigned int i = 0; i < 10; ++i) m[i] = c; }

        SymetricMatrix(float m11, float m12, float m13, float m14,
            float m22, float m23, float m24,
            float m33, float m34,
            float m44) {
			m[0] = m11;  m[1] = m12;  m[2] = m13;  m[3] = m14;
			m[4] = m22;  m[5] = m23;  m[6] = m24;
			m[7] = m33;  m[8] = m34;
			m[9] = m44;
		}

        SymetricMatrix(float a, float b, float c, float d) {
			m[0] = a * a;  m[1] = a * b;  m[2] = a * c;  m[3] = a * d;
			m[4] = b * b;  m[5] = b * c;  m[6] = b * d;
			m[7] = c * c;  m[8] = c * d;
			m[9] = d * d;
		}

        float operator[](int c) const { return m[c]; }

		// Determinant
        float det(int a11, int a12, int a13,
			int a21, int a22, int a23,
			int a31, int a32, int a33) {
            float det = m[a11] * m[a22] * m[a33] + m[a13] * m[a21] * m[a32] + m[a12] * m[a23] * m[a31]
				- m[a13] * m[a22] * m[a31] - m[a11] * m[a23] * m[a32] - m[a12] * m[a21] * m[a33];
			return det;
		}

		const SymetricMatrix operator+(const SymetricMatrix& n) const {
			return SymetricMatrix(m[0] + n[0], m[1] + n[1], m[2] + n[2], m[3] + n[3],
				m[4] + n[4], m[5] + n[5], m[6] + n[6],
				m[7] + n[7], m[8] + n[8],
				m[9] + n[9]);
		}

		SymetricMatrix& operator+=(const SymetricMatrix& n) {
			m[0] += n[0];   m[1] += n[1];   m[2] += n[2];   m[3] += n[3];
			m[4] += n[4];   m[5] += n[5];   m[6] += n[6];   m[7] += n[7];
			m[8] += n[8];   m[9] += n[9];
			return *this;
		}

        float m[10];
	};

	class Mesh_Simplifier {
	public:
        struct Triangle { int v[3]; float err[4]; int deleted, dirty; Vector3 n; };
		struct Vertex { Vector3 p; int tstart, tcount; SymetricMatrix q; int border; };
		struct Ref { int tid, tvertex; };
		std::vector<Triangle> triangles;
		std::vector<Vertex> vertices;
		std::vector<Ref> refs;

		// target_count  : target nr. of triangles
		// agressiveness : sharpness to increase the threshold.
		//                 5..8 are good numbers
		//                 more iterations yield higher quality
		//
        void simplify_mesh(int target_count, float agressiveness = 7)
		{
			for (unsigned int i = 0; i < triangles.size(); ++i) triangles[i].deleted = 0;

			int deleted_triangles = 0;
			std::vector<int> deleted0, deleted1;
			int triangle_count = triangles.size();

			for (int iteration = 0; iteration < 100; ++iteration)
			{
				if (triangle_count - deleted_triangles <= target_count)break;

				// update mesh once in a while
				if (iteration % 5 == 0)
				{
					update_mesh(iteration);
				}

				// clear dirty flag
				for (unsigned int i = 0; i < triangles.size(); ++i) triangles[i].dirty = 0;

				//
				// All triangles with edges below the threshold will be removed
				//
				// The following numbers works well for most models.
				// If it does not, try to adjust the 3 parameters
				//
                float threshold = 0.000000001*pow(float(iteration + 3), agressiveness);

				// remove vertices & mark deleted triangles			
				for (unsigned int i = 0; i < triangles.size(); ++i)
				{
					Triangle &t = triangles[i];
					if (t.err[3] > threshold) continue;
					if (t.deleted) continue;
					if (t.dirty) continue;

					for (int j = 0; j < 3; j++)
						if (t.err[j] < threshold)
						{
							int i0 = t.v[j]; Vertex &v0 = vertices[i0];
							int i1 = t.v[(j + 1) % 3]; Vertex &v1 = vertices[i1];

							// Border check
							if (v0.border != v1.border)  continue;

							// Compute vertex to collapse to
							Vector3 p;
							calculate_error(i0, i1, p);

							deleted0.resize(v0.tcount); // normals temporarily
							deleted1.resize(v1.tcount); // normals temporarily

							// don't remove if flipped
							if (flipped(p, i0, i1, v0, v1, deleted0)) continue;
							if (flipped(p, i1, i0, v1, v0, deleted1)) continue;

							// not flipped, so remove edge												
							v0.p = p;
							v0.q = v1.q + v0.q;
							int tstart = refs.size();

							update_triangles(i0, v0, deleted0, deleted_triangles);
							update_triangles(i0, v1, deleted1, deleted_triangles);

							int tcount = refs.size() - tstart;

							if (tcount <= v0.tcount)
							{
								// save ram
								if (tcount)memcpy(&refs[v0.tstart], &refs[tstart], tcount * sizeof(Ref));
							}
							else
								// append
								v0.tstart = tstart;

							v0.tcount = tcount;
							break;
						}
					// done?
					if (triangle_count - deleted_triangles <= target_count)break;
				}
			}
			compact_mesh();
		}

		// Check if a triangle flips when this edge is removed
		bool flipped(Vector3 p, int i0, int i1, Vertex &v0, Vertex &v1, std::vector<int> &deleted)
		{
			int bordercount = 0;
			for (unsigned int k = 0; k < v0.tcount; ++k)
			{
				Triangle &t = triangles[refs[v0.tstart + k].tid];
				if (t.deleted)continue;

				int s = refs[v0.tstart + k].tvertex;
				int id1 = t.v[(s + 1) % 3];
				int id2 = t.v[(s + 2) % 3];

				if (id1 == i1 || id2 == i1) // delete ?
				{
					bordercount++;
					deleted[k] = 1;
					continue;
				}
				Vector3 d1 = vertices[id1].p - p; d1.normalize();
				Vector3 d2 = vertices[id2].p - p; d2.normalize();
				if (fabs(d1.dot(d2)) > 0.999) return true;
				Vector3 n = d1.cross(d2);
				n.normalize();
				deleted[k] = 0;
				if (n.dot(t.n) < 0.2) return true;
			}
			return false;
		}

		// Update triangle connections and edge error after a edge is collapsed
		void update_triangles(int i0, Vertex &v, std::vector<int> &deleted, int &deleted_triangles)
		{
			Vector3 p;
			for (unsigned int k = 0; k < v.tcount; ++k)
			{
				Ref &r = refs[v.tstart + k];
				Triangle &t = triangles[r.tid];
				if (t.deleted)continue;
				if (deleted[k])
				{
					t.deleted = 1;
					deleted_triangles++;
					continue;
				}
				t.v[r.tvertex] = i0;
				t.dirty = 1;
				t.err[0] = calculate_error(t.v[0], t.v[1], p);
				t.err[1] = calculate_error(t.v[1], t.v[2], p);
				t.err[2] = calculate_error(t.v[2], t.v[0], p);
				t.err[3] = std::fmin(t.err[0], std::fmin(t.err[1], t.err[2]));
				refs.push_back(r);
			}
		}

		// compact triangles, compute edge error and build reference list
		void update_mesh(int iteration) {
			if (iteration > 0) // compact triangles
			{
				int dst = 0;
				for (unsigned int i = 0; i < triangles.size(); ++i)
					if (!triangles[i].deleted)
					{
						triangles[dst++] = triangles[i];
					}
				triangles.resize(dst);
			}
			//
			// Init Quadrics by Plane & Edge Errors
			//
			// required at the beginning ( iteration == 0 )
			// recomputing during the simplification is not required,
			// but mostly improves the result for closed meshes
			//
			if (iteration == 0) {
				for (unsigned int i = 0; i < vertices.size(); ++i)
					vertices[i].q = SymetricMatrix(0.0);

				for (unsigned int i = 0; i < triangles.size(); ++i) {
					Triangle &t = triangles[i];
					Vector3 n, p[3];
					for (int j = 0; j < 3; ++j) p[j] = vertices[t.v[j]].p;
					n = (p[1] - p[0]).cross(p[2] - p[0]);
					n.normalize();
					t.n = n;
					for (int j = 0; j < 3; ++j) vertices[t.v[j]].q =
						vertices[t.v[j]].q + SymetricMatrix(n.x, n.y, n.z, -n.dot(p[0]));
				}
				for (unsigned int i = 0; i < triangles.size(); ++i) {
					// Calc Edge Error
					Triangle &t = triangles[i]; Vector3 p;
					for (int j = 0; j < 3; ++j) t.err[j] = calculate_error(t.v[j], t.v[(j + 1) % 3], p);
					t.err[3] = std::fmin(t.err[0], std::fmin(t.err[1], t.err[2]));
				}
			}

			// Init Reference ID list	
			for (unsigned int i = 0; i < vertices.size(); ++i)
			{
				vertices[i].tstart = 0;
				vertices[i].tcount = 0;
			}
			for (unsigned int i = 0; i < triangles.size(); ++i)
			{
				Triangle &t = triangles[i];
				for (int j = 0; j < 3; ++j) vertices[t.v[j]].tcount++;
			}
			int tstart = 0;
			for (unsigned int i = 0; i < vertices.size(); ++i)
			{
				Vertex &v = vertices[i];
				v.tstart = tstart;
				tstart += v.tcount;
				v.tcount = 0;
			}

			// Write References
			refs.resize(triangles.size() * 3);
			for (unsigned int i = 0; i < triangles.size(); ++i)
			{
				Triangle &t = triangles[i];
				for (int j = 0; j < 3; ++j)
				{
					Vertex &v = vertices[t.v[j]];
					refs[v.tstart + v.tcount].tid = i;
					refs[v.tstart + v.tcount].tvertex = j;
					v.tcount++;
				}
			}

			// Identify boundary : vertices[].border=0,1 
			if (iteration == 0)
			{
				std::vector<int> vcount, vids;

				for (unsigned int i = 0; i < vertices.size(); ++i)
					vertices[i].border = 0;

				for (unsigned int i = 0; i < vertices.size(); ++i)
				{
					Vertex &v = vertices[i];
					vcount.clear();
					vids.clear();
					for (int j = 0; j < v.tcount; j++)
					{
						int k = refs[v.tstart + j].tid;
						Triangle &t = triangles[k];
						for (int k = 0; k < 3; ++k)
						{
							int ofs = 0, id = t.v[k];
							while (ofs < vcount.size())
							{
								if (vids[ofs] == id)break;
								ofs++;
							}
							if (ofs == vcount.size())
							{
								vcount.push_back(1);
								vids.push_back(id);
							}
							else
								vcount[ofs]++;
						}
					}
					for (int j = 0; j < vcount.size(); ++j)
						if (vcount[j] == 1)
							vertices[vids[j]].border = 1;
				}
			}
		}

		// Finally compact mesh before exiting
		void compact_mesh() {
			int dst = 0;
			for (unsigned int i = 0; i < vertices.size(); ++i)
			{
				vertices[i].tcount = 0;
			}
			for (unsigned int i = 0; i < triangles.size(); ++i)
				if (!triangles[i].deleted) {
					Triangle &t = triangles[i];
					triangles[dst++] = t;
					for (int j = 0; j < 3; ++j) vertices[t.v[j]].tcount = 1;
				}
			triangles.resize(dst);
			dst = 0;
			for (unsigned int i = 0; i < vertices.size(); ++i)
				if (vertices[i].tcount) {
					vertices[i].tstart = dst;
					vertices[dst].p = vertices[i].p;
					dst++;
				}
			for (unsigned int i = 0; i < triangles.size(); ++i) {
				Triangle &t = triangles[i];
				for (int j = 0; j < 3; ++j) t.v[j] = vertices[t.v[j]].tstart;
			}
			vertices.resize(dst);
		}

		// Error between vertex and Quadric
        float vertex_error(SymetricMatrix q, float x, float y, float z) {
			return   q[0] * x*x + 2 * q[1] * x*y + 2 * q[2] * x*z + 2 * q[3] * x + q[4] * y*y
				+ 2 * q[5] * y*z + 2 * q[6] * y + q[7] * z*z + 2 * q[8] * z + q[9];
		}

		// Error for one edge
        float calculate_error(int id_v1, int id_v2, Vector3 &p_result) {
			// compute interpolated vertex 
			SymetricMatrix q = vertices[id_v1].q + vertices[id_v2].q;
			bool   border = vertices[id_v1].border & vertices[id_v2].border;
            float error = 0;
            float det = q.det(0, 1, 2, 1, 4, 5, 2, 5, 7);

			if (det != 0 && !border) {
				// q_delta is invertible
				p_result.x = -1 / det * (q.det(1, 2, 3, 4, 5, 6, 5, 7, 8));	// vx = A41/det(q_delta) 
				p_result.y = 1 / det * (q.det(0, 2, 3, 1, 5, 6, 2, 7, 8));	// vy = A42/det(q_delta) 
				p_result.z = -1 / det * (q.det(0, 1, 3, 1, 4, 6, 2, 5, 8));	// vz = A43/det(q_delta) 
				error = vertex_error(q, p_result.x, p_result.y, p_result.z);
			}
			else {
				// det = 0 -> try to find best result
				Vector3 p1 = vertices[id_v1].p;
				Vector3 p2 = vertices[id_v2].p;
				Vector3 p3 = (p1 + p2) / 2;
                float error1 = vertex_error(q, p1.x, p1.y, p1.z);
                float error2 = vertex_error(q, p2.x, p2.y, p2.z);
                float error3 = vertex_error(q, p3.x, p3.y, p3.z);
				error = std::fmin(error1, std::fmin(error2, error3));
				if (error1 == error) p_result = p1;
				if (error2 == error) p_result = p2;
				if (error3 == error) p_result = p3;
			}
			return error;
		}
	};
};
