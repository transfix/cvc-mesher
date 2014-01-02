TEMPLATE = lib
CONFIG  += qt warn_on staticlib create_prl
#TARGET  += FastContouring 

INCLUDEPATH += ../VolMagick

QMAKE_CXXFLAGS += $$(CPPFLAGS)
QMAKE_LFLAGS += $$(LDFLAGS)

# Input
HEADERS = \
	ContourGeometry.h \
	cubes.h \
	FastContouring.h \
	MarchingCubesBuffers.h \
	Matrix.h \
	Quaternion.h \
	Ray.h \
	Tuple.h \
	Vector.h

SOURCES = \
	ContourGeometry.cpp \
	FastContouring.cpp \
	MarchingCubesBuffers.cpp \
	Matrix.cpp \
	Quaternion.cpp \
	Ray.cpp  \
	Tuple.cpp \
	Vector.cpp
