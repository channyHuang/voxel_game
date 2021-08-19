# voxel_game

using voxel to edite game terrain

input: sdf value stored by VoxelBuffer

output: .obj mesh file

## Overview

Using class VoxelBuffer in voxel tool to store data, a plugin of godot engine, but simplify it because this project is mainly about sdf data generation and mesh generation. So this project do NOT consider memory cases.

使用了godot引擎插件voxel tool中的VoxelBuffer类来存储数据，但做了简化。因为本项目主要研究sdf数据的生成和网格生成，不考虑内存相关逻辑。

## Mesh Generation

### Marching Cube

### Dual Marching Cube

### Transvoxel

### Naive Surface Nets

every voxel generate less than or equal to one vertex

### Dual Contouring

need normals before generate surface

## Reference

[godot engine](https://godotengine.org/)  (MIT license) 

[voxel tool](https://voxel-tools.readthedocs.io/en/latest/) (GNU license) godot plugin, using its VoxelBuffer class

[naive surface nets](https://github.com/Q-Minh/naive-surface-nets)

[dual contouring](https://github.com/emilk/Dual-Contouring)

[fast noise lite](https://github.com/Auburn/FastNoiseLite/wiki/Documentation)