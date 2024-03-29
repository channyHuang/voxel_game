#include($$PWD/dualContouring/dualContouring.pri)
include($$PWD/naiveSurfaceNets/naiveSurfaceNets.pri)

HEADER_FILES = $$files($$PWD/*.h, false)
HEADER_SOURCE_FILES = $$files($$PWD/*.hpp, false)
SOURCE_FILES = $$files($$PWD/*.cpp, false)

INCLUDEPATH += $$PWD/../

HEADERS += $${HEADER_FILES} \
    $${HEADER_SOURCE_FILES}

SOURCES += $${SOURCE_FILES}
