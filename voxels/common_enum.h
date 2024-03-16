#ifndef COMMON_ENUM_H
#define COMMON_ENUM_H

namespace TerrainEnum {

enum Notification_Event {
    Notification_Enter = 1,
    Notification_Exit,
    Notification_Process
};

enum MaterialType {
    AIR = 0,
    BASALT,
    BRICK,
    GRASS,
    GRASSLAND,
    ICE,
    LAVA,
    MUD,
    ROCK,
    SAND,
    SANDSTONE,
    SNOW,
    SOIL,
    STONE,
    WATER,
    WOOD,
    ROAD,
    Material_MAX
};

}

#endif // COMMON_ENUM_H
