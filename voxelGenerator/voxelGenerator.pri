HEADER_FILES = $$files($$PWD/*.h, false)
SOURCE_FILES = $$files($$PWD/*.cpp, false)

HEADERS += $${HEADER_FILES}
SOURCES += $${SOURCE_FILES}
