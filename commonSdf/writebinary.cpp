#include "voxel_buffer.h"
#include "sdf3.h"

#include <fstream>

bool readBinvox(std::string sfilename, VoxelBuffer &voxels) {
    std::ifstream ifs(sfilename);
    if (!ifs.is_open()) {
        return false;
    }
    Vector3i vdim;
    char c[3];
    return true;
    ifs.read(c, 3);
    vdim.x = c[0];
    vdim.y = c[1];
    vdim.z = c[2];

    char csdf;
    for (int i = -16; i < 15; ++i) {
        for (int j = -16; j < 15; ++j) {
            for (int k = -16; k < 15; ++k) {
                ifs.read(&csdf, 1);
                uint8_t sdf = static_cast<uint8_t>(csdf);
                Vector3i vposi = Vector3i(i, j, k);
                voxels.set_voxel(sdf, vposi);
            }
        }
    }
    ifs.close();
}

bool writeBinvox(std::string sfilename, const VoxelBuffer &voxels) {
    std::ofstream ofs(sfilename, std::ios::out | std::ios::binary);
    Vector3i vdim = voxels.get_size();
    ofs.write(reinterpret_cast<char*>(&vdim.x), 1);
    ofs.write(reinterpret_cast<char*>(&vdim.y), 1);
    ofs.write(reinterpret_cast<char*>(&vdim.z), 1);

    for (int i = 0; i < vdim.x; ++i) {
        for (int j = 0; j < vdim.y; ++j) {
            for (int k = 0; k < vdim.z; ++k) {
                uint8_t sdfint = static_cast<uint8_t>(voxels.get_voxel(i, j, k));
                ofs.write(reinterpret_cast<char*>(&sdfint), 1);
            }
        }
    }
    ofs.close();
    return true;
}

void generate() {
    Sdf3 sdf3;
    for (uint8_t i = 0; i < static_cast<uint8_t>(Sdf3::SdfType::Sdf_Count); i++) {
        Sdf3::SdfType type = static_cast<Sdf3::SdfType>(i);
        switch (type) {
        case Sdf3::SdfType::Sdf_Sphere:
            sdf3.Sphere(10, Vector3(0, 0, 0));
            break;
        case Sdf3::SdfType::Sdf_Slab:
            sdf3.Slab(Vector3(0, -10, 0), Vector3(0, 10, 0), 0.0f);
            break;
        case Sdf3::SdfType::Sdf_Box:
            sdf3.Box(Vector3(10, 10, 10), Vector3(0, 0, 0));
            break;
        case Sdf3::SdfType::Sdf_RoundBox:
            sdf3.RoundBox(Vector3(10, 10, 10), 3.f);
            break;
        case Sdf3::SdfType::Sdf_WireframeBox:
            sdf3.WireframeBox(Vector3(10, 10, 10), 3.f);
            break;
        case Sdf3::SdfType::Sdf_Torus:
            sdf3.Torus(3.f, 5.f);
            break;
        case Sdf3::SdfType::Sdf_Capsule:
            sdf3.Capsule(Vector3(0, -10, 0), Vector3(0, 10, 0), 3.f);
            break;
        case Sdf3::SdfType::Sdf_Cylinder:
            sdf3.Cylinder(3.f);
            break;
        case Sdf3::SdfType::Sdf_CappedCone:
            sdf3.CappedCone(Vector3(0, -10, 0), Vector3(0, 10, 0), 3.f, 5.f);
            break;
        case Sdf3::SdfType::Sdf_RoundedCone:
            sdf3.RoundedCone(3.f, 5.f, 10.f);
            break;
        case Sdf3::SdfType::Sdf_Ellipsoid:
            sdf3.Ellipsoid(10.f);
            break;
        case Sdf3::SdfType::Sdf_Pyramid:
            sdf3.Pyramid(10.f);
            break;
        case Sdf3::SdfType::Sdf_Tetrahedron:
            sdf3.Tetrahedron(10);
            break;
            break;
        case Sdf3::SdfType::Sdf_Octahedron:
            sdf3.Octahedron(10.f);
            break;
        case Sdf3::SdfType::Sdf_Dodecahedron:
            sdf3.Dodecahedron(10.f);
            break;
        case Sdf3::SdfType::Sdf_TriPrism:
            sdf3.TriPrism(Vector2(10.f, 10.f));
            break;
        default:
            break;
        }

        std::ofstream ofs(Sdf3::vSdfTypeString[i] + ".voxel", std::ios::out | std::ios::binary);
        Vector3i vdim = Vector3i(32, 32, 32);
        ofs.write(reinterpret_cast<char*>(&vdim.x), 1);
        ofs.write(reinterpret_cast<char*>(&vdim.y), 1);
        ofs.write(reinterpret_cast<char*>(&vdim.z), 1);

        for (int i = -16; i < 15; ++i) {
            for (int j = -16; j < 15; ++j) {
                for (int k = -16; k < 15; ++k) {
                    Vector3 vposi = Vector3(i, j, k);
                    real sdf = sdf3.getFun()(vposi);
                    uint8_t sdfint = Math::clamp(static_cast<int>(128.f * sdf + 128.f), 0, 0xfe);

                    ofs.write(reinterpret_cast<char*>(&sdfint), 1);
                }
            }
        }
        ofs.close();
    }
    if (false) {
        std::string inputstr = "";
        std::ifstream ifs("voxelbinary/" + inputstr + ".voxel", std::ios::in | std::ios::binary);
        Vector3i vdim;
        char c[3];
        ifs.read(c, 3);
        vdim.x = c[0];
        vdim.y = c[1];
        vdim.z = c[2];

        char csdf;
        for (int i = -16; i < 15; ++i) {
            for (int j = -16; j < 15; ++j) {
                for (int k = -16; k < 15; ++k) {
                    ifs.read(&csdf, 1);
                    Real sdf = (static_cast<uint8_t>(csdf) - 0x7f) * 1.f / 0x7f;
                    Vector3 vposi = Vector3(i, j, k);
                }
            }
        }
        ifs.close();
    }
}
