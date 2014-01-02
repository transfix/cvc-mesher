TEMPLATE = lib
#TEMPLATE = app
CONFIG  += warn_on staticlib rtti exceptions
CONFIG -= qt
#TARGET  += VolMagick 

QMAKE_CXXFLAGS += $$(CPPFLAGS)
QMAKE_LFLAGS += $$(LDFLAGS)

SOURCES = VolMagick.cpp \
	RawIV.cpp \
	RawV.cpp \ 
	MRC.cpp \
	INR.cpp \
	Spider.cpp \
	VolumeCache.cpp \
	AnisotropicDiffusion.cpp \
	BilateralFilter.cpp \
	ContrastEnhancement.cpp #main.cpp

HEADERS = BoundingBox.h \
	Dimension.h \
	Exceptions.h \
	VolMagick.h \
	endians.h \
	VolumeCache.h

solaris-cc | solaris-cc-64 | solaris-g++ | solaris-g++-64 {
DEFINES += __SOLARIS__
} 

win32-g++ | win32-msvc {
DEFINES += __WINDOWS__
}
