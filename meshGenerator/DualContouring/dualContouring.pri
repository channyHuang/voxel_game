HEADER_FILES = $$files($$PWD/*.h, false)
HEADER_SOURCE_FILES = $$files($$PWD/*.hpp, false)
SOURCE_FILES = $$files($$PWD/*.cpp, false)

HEADERS += $${HEADER_FILES} \
    $${HEADER_SOURCE_FILES}
SOURCES += $${SOURCE_FILES}
