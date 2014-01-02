# General settings

TEMPLATE = lib
CONFIG += create_prl warn_off staticlib
#TARGET  += contour
DEFINES += _BOOL  SP2 USEDICT NDEBUG

# Don't define _LITTLE_ENDIAN on MacOS X (it's big endian)
!macx {
	DEFINES += _LITTLE_ENDIAN
} else {
	DEFINES += MACOS_X
}

win32-g++ {
	QMAKE_CXXFLAGS += -fpermissive
}

# Input

SOURCES =  \
		bin.cpp \
		contour.cpp \
		inttree.cpp \
		seedall.cpp \
		bucketsearch.cpp \
		cubes.cpp \
		ipqueue.cpp \
		seedcells.cpp \
		cellqueue.cpp \
		data.cpp \
		iqueue.cpp \
		seedchkr2.cpp \
		conplot2d.cpp \
		datareg2.cpp \
		seedchkr3.cpp \
		conplot3d.cpp \
		datareg3.cpp \
		queue.cpp \
		seeddirreg3.cpp \
		conplot.cpp \
		dataslc.cpp \
		range.cpp \
		segtree.cpp \
		conplot_p.cpp \
		datavol.cpp \
		rangeprop.cpp \
		shelf.cpp \
		conplotreg2.cpp \
		dirseeds.cpp \
		rangesweep.cpp \
		spectrumtest.cpp \
		conplotreg3.cpp \
		dirseedsreg2.cpp \
		regprop2.cpp \
		squeue.cpp \
		contour2d.cpp \
		edgehash.cpp \
		regprop.cpp \
		contour3d.cpp \
		hash.cpp \
		respprop2.cpp \
		dict.c

HEADERS =  \
		basic.h \
		contour3d.h \
		dirseeds.h \
		regprop.h \
		bin.h \
		contour.h \
		dirseedsreg2.h \
		respprop2.h \
		bucketsearch.h \
		cubes.h \
		edgehash.h \
		seedall.h \
		cellqueue.h \
		data.h \
		endian_io.h \
		seedcells.h \
		cellsearch.h \
		datareg2.h \
		hash.h \
		seedchkr2.h \
		commdata.h \
		datareg3.h \
		inttree.h \
		seedchkr3.h \
		compute.h \
		dataset.h \
		ipqueue.h \
		seeddirreg3.h \
		conplot2d.h \
		datasetreg2.h \
		iqueue.h \
		segtree.h \
		conplot3d.h \
		datasetreg3.h \
		queue.h \
		shelf.h \
		conplot.h \
		datasetslc.h \
		range.h \
		squeue.h \
		conplotreg2.h \
		datasetvol.h \
		rangeprop.h \
		utilities.h \
		conplotreg3.h \
		dataslc.h \
		rangesweep.h \
		vtkmarchingcubescases.h \
		contour2d.h \
		datavol.h \
		regprop2.h \
		dict.h 

