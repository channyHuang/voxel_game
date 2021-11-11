HEADER_FILES = $$files($$PWD/*.h, false)
SOURCE_FILES = $$files($$PWD/*.cpp, false)

HEADERS += $${HEADER_FILES} \
    $$PWD/terraingenerator_roblox.h
SOURCES += $${SOURCE_FILES} \
    $$PWD/terraingenerator_roblox.cpp
