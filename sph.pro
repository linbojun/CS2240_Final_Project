QT += core gui opengl widgets

TARGET = simulation
TEMPLATE = app

QMAKE_CXXFLAGS += -mstackrealign
CONFIG += c++17

unix:!macx {
    LIBS += -lGLU
}
win32 {
    DEFINES += GLEW_STATIC
    LIBS += -lopengl32 -lglu32
}

SOURCES += \
    libs/glew-1.10.0/src/glew.c \
    src/SPH.cpp \
    src/graphics/Shader.cpp \
    src/graphics/GraphicsDebug.cpp \
    src/graphics/shape.cpp \
    src/graphics/camera.cpp \
    src/graphics/MeshLoader.cpp \
    src/magneticSPH.cpp \
    src/main.cpp \
    src/mainwindow.cpp \
    src/simdump.cpp \
    src/simulation.cpp \
    src/view.cpp \
    src/viewformat.cpp

HEADERS += \
    src/SPH.h \
    src/graphics/Shader.h \
    src/graphics/ShaderAttribLocations.h \
    src/graphics/GraphicsDebug.h \
    src/graphics/shape.h \
    src/graphics/camera.h \
    src/mainwindow.h \
    src/magneticSPH.h \
    src/simdump.h \
    src/simulation.h \
    src/view.h \
    src/viewformat.h \
    ui_mainwindow.h \
    src/graphics/MeshLoader.h

FORMS += src/mainwindow.ui

RESOURCES += \
    res/shaders/shaders.qrc

DISTFILES += \
    res/shaders/shader.vert \
    res/shaders/shader.frag

INCLUDEPATH += src libs glm libs/glew-1.10.0/include libs/Eigen/
DEPENDPATH += src libs glm libs/glew-1.10.0/include libs/Eigen/

