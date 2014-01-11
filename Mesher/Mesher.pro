TEMPLATE = app
CONFIG += warn_on thread rtti exceptions console
CONFIG -= qt
TARGET = mesher

SOURCES = main.cpp mesher.cpp
HEADERS = mesher.h

QMAKE_CXXFLAGS += $$(CPPFLAGS)
QMAKE_LFLAGS += $$(LDFLAGS)

INCLUDEPATH += ../LBIE ../contour ../VolMagick ../FastContouring
LIBS += ../LBIE/libLBIE.a \
	../VolMagick/libVolMagick.a \
	../contour/libcontour.a \
	../FastContouring/libFastContouring.a \
	-lboost_program_options

