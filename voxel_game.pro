FQT += core gui widgets

CONFIG += c++2a
CONFIG -= app_bundle
CONFIG += sdk_no_version_check

# The following define makes your compiler emit warnings if you use
# any Qt feature that has been marked deprecated (the exact warnings
# depend on your compiler). Please consult the documentatFion of the
# deprecated API in order to know how to port your code away from it.
DEFINES += QT_DEPRECATED_WARNINGS

# You can also make your code fail to compile if it uses deprecated APIs.
# In order to do so, uncomment the following line.
# You can also select to disable deprecated APIs only up to a certain version of Qt.
#DEFINES += QT_DISABLE_DEPRECATED_BEFORE=0x060000    # disables all the APIs deprecated before Qt 6.0.0

INCLUDEPATH += $$PWD/../commonClass/
include($$PWD/../commonClass/commonMath/commonMath.pri)
include($$PWD/../commonClass/commonGeometry/commonGeometry.pri)
include($$PWD/../commonClass/commonFunc/commonFunc.pri)

include($$PWD/commonSdf/commonSdf.pri)
include($$PWD/biomes/biomes.pri)

include($$PWD/voxels/voxels.pri)
include($$PWD/meshGenerator/meshGenerator.pri)
include($$PWD/voxelGenerator/voxelGenerator.pri)
include($$PWD/terrains/terrains.pri)

include($$PWD/renderControl/renderControl.pri)

INCLUDEPATH += $$PWD

HEADER_FILES = $$files($$PWD/*.h, false)
SOURCE_FILES = $$files($$PWD/*.cpp, false)
RESOURCE_FILES = $$files($$PWD/*.qrc, false)

HEADERS += $${HEADER_FILES}
SOURCES += $${SOURCE_FILES}
RESOURCES += $${RESOURCE_FILES}


# Default rules for deployment.
qnx: target.path = /tmp/$${TARGET}/bin
else: unix:!android: target.path = /opt/$${TARGET}/bin
!isEmpty(target.path): INSTALLS += target
