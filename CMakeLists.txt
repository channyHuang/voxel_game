cmake_minimum_required(VERSION 3.14)

project(3DTerrain LANGUAGES CXX)

set(CMAKE_AUTOUIC ON)
set(CMAKE_AUTOMOC ON)
set(CMAKE_AUTORCC ON)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(QT NAMES Qt6 Qt5 REQUIRED COMPONENTS Core)
find_package(Qt${QT_VERSION_MAJOR} REQUIRED COMPONENTS Core)

set(osg_path /Users/channyhuang/Documents/thirdlib/OpenSceneGraph/include)
set(osg_lib /Users/channyhuang/Documents/thirdlib/OpenSceneGraph/build/lib)

include_directories(${osg_path})
include_directories("/usr/local/include")

set(IMGUI_PATH ./../commonClass/imgui)
include_directories("./../commonClass/imgui")

link_directories("/usr/local/lib")
link_directories("/usr/local/lib")
link_directories("/System/Library/Frameworks")
link_directories(${osg_lib})

add_definitions(-DGLFW_STATIC)

file(GLOB ALL_SOURCE ./*.cpp)
file(GLOB ALL_HEAD ./*.h)
file(GLOB IMGUI_FILE ${IMGUI_PATH}/*.*)


file(GLOB VOXEL_GEN_FILE ./voxelGenerator/*)
file(GLOB VOXEL_FILE ./voxels/*)
file(GLOB TERRAIN_FILE ./terrains/*)
file(GLOB MESH_FILE ./meshGenerator/*)
file(GLOB MESH_GEN_FILE ./meshGenerator/naiveSurfaceNets/*)

file(GLOB COMMON_FILE ./../commonClass/*/*)

include_directories("./../commonClass")
include_directories("./")

add_executable(3DTerrain
    ${ALL_SOURCE}
    ${ALL_HEAD}
    ${IMGUI_FILE}

    ${VOXEL_GEN_FILE}
    ${VOXEL_FILE}
    ${TERRAIN_FILE}
    ${MESH_FILE}
    ${MESH_GEN_FILE}

    ${COMMON_FILE}
)

target_link_libraries(3DTerrain #[[Qt${QT_VERSION_MAJOR}::Core]]
    libglfw3.a
    libGLEW.a
    "-framework Cocoa"
    "-framework OpenGL"
    "-framework CoreVideo"
    "-framework IOKit"
    osg
    osgViewer
    osgUtil
    osgGA
    OpenThreads

)

install(TARGETS 3DTerrain
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})
