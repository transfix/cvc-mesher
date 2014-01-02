TEMPLATE = lib
CONFIG  += qt warn_on staticlib create_prl
#TARGET  += LBIE 

INCLUDEPATH += ../VolMagick ../contour ../FastContouring

QMAKE_CXXFLAGS += $$(CPPFLAGS)
QMAKE_LFLAGS += $$(LDFLAGS)

# Input
HEADERS = \
		LBIE_Mesher.h \
                e_face.h \
                LBIE_geoframe.h \
                normalspline.h \
                octree.h \
                pcio.h

SOURCES = \
		LBIE_Mesher.cpp \
		e_face.cpp \
		LBIE_geoframe.cpp \
		hexa.cpp \
		normalspline.cpp \
		octree.cpp \
		pcio.cpp \
		tetra.cpp
