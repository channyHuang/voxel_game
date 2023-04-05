QT       += core openglwidgets

greaterThan(QT_MAJOR_VERSION, 4): QT += widgets

CONFIG += c++11

# The following define makes your compiler emit warnings if you use
# any Qt feature that has been marked deprecated (the exact warnings
# depend on your compiler). Please consult the documentation of the
# deprecated API in order to know how to port your code away from it.
DEFINES += QT_DEPRECATED_WARNINGS

include($$PWD/../renderCommon/renderCommon.pri)
include($$PWD/../uvUnwrap/uvUnwrap.pri)
INCLUDEPATH += $$PWD/../renderCommon


SOURCE_FILES = $$files($${PWD}/*.cpp, false)
HEADER_FILES = $$files($${PWD}/*.h, false)
RESOURCE_FILES = $$files($${PWD}/*.qrc, false)

SOURCES += $${SOURCE_FILES}

HEADERS += $${HEADER_FILES}

DEFINES += PRO_PATH=$${PWD}

RESOURCES += $${RESOURCE_FILES}

# Default rules for deployment.
qnx: target.path = /tmp/$${TARGET}/bin
else: unix:!android: target.path = /opt/$${TARGET}/bin
!isEmpty(target.path): INSTALLS += target
