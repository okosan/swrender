

INCLUDEPATH += $$PWD
SOURCES += $$PWD/rasterwindow.cpp
HEADERS += $$PWD/rasterwindow.h

SOURCES += \
    main.cpp

target.path = $$PWD/render
INSTALLS += target

CONFIG += c++11
