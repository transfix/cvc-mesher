# Microsoft Developer Studio Project File - Name="contour" - Package Owner=<4>
# Microsoft Developer Studio Generated Build File, Format Version 6.00
# ** DO NOT EDIT **

# TARGTYPE "Win32 (x86) Static Library" 0x0104

CFG=contour - Win32 Release
!MESSAGE This is not a valid makefile. To build this project using NMAKE,
!MESSAGE use the Export Makefile command and run
!MESSAGE 
!MESSAGE NMAKE /f "contour.mak".
!MESSAGE 
!MESSAGE You can specify a configuration when running NMAKE
!MESSAGE by defining the macro CFG on the command line. For example:
!MESSAGE 
!MESSAGE NMAKE /f "contour.mak" CFG="contour - Win32 Release"
!MESSAGE 
!MESSAGE Possible choices for configuration are:
!MESSAGE 
!MESSAGE "contour - Win32 Release" (based on "Win32 (x86) Static Library")
!MESSAGE "contour - Win32 Debug" (based on "Win32 (x86) Static Library")
!MESSAGE 

# Begin Project
# PROP AllowPerConfigDependencies 0
# PROP Scc_ProjName ""
# PROP Scc_LocalPath ""
CPP=cl.exe
RSC=rc.exe

!IF  "$(CFG)" == "contour - Win32 Release"

# PROP BASE Use_MFC 0
# PROP BASE Use_Debug_Libraries 0
# PROP BASE Output_Dir "Release"
# PROP BASE Intermediate_Dir "Release"
# PROP BASE Target_Dir ""
# PROP Use_MFC 0
# PROP Use_Debug_Libraries 0
# PROP Output_Dir "Release"
# PROP Intermediate_Dir "Release"
# PROP Target_Dir ""
# ADD CPP /nologo /W3 /GX /O1 /I "$(QTDIR)\include" /I "D:\users\transfix\VolumeRover\libcontour" /I "d:\Qt\3.3.3Educational\mkspecs\win32-msvc" /D "WIN32" /D "NDEBUG" /D "_LIB" /D "UNICODE" /D "_BOOL" /D "SP2" /D "USEDICT" /D "_LITTLE_ENDIAN" /D "QT_DLL" /D "QT_THREAD_SUPPORT" /D "QT_NO_DEBUG" /FD -Zm200 /c
# ADD BASE RSC /l 0x409
# ADD RSC /l 0x409 /d "NDEBUG"
BSC32=bscmake.exe
# ADD BSC32 /nologo
LIB32=link.exe -lib
# ADD LIB32 /nologo

!ELSEIF  "$(CFG)" == "contour - Win32 Debug"

# PROP BASE Use_MFC 0
# PROP BASE Use_Debug_Libraries 1
# PROP BASE Output_Dir "Debug"
# PROP BASE Intermediate_Dir "Debug"
# PROP BASE Target_Dir ""
# PROP Use_MFC 0
# PROP Use_Debug_Libraries 1
# PROP Output_Dir "Debug"
# PROP Intermediate_Dir "Debug"
# PROP Target_Dir ""
# ADD CPP /nologo /MDd /W3 /GX /Zi /Od /I "$(QTDIR)\include" /I "D:\users\transfix\VolumeRover\libcontour" /I "d:\Qt\3.3.3Educational\mkspecs\win32-msvc" /D "WIN32" /D "_DEBUG" /D "_LIB" /D "UNICODE" /D "_BOOL" /D "SP2" /D "USEDICT" /D "NDEBUG" /D "_LITTLE_ENDIAN" /D "QT_DLL" /D "QT_THREAD_SUPPORT" /FD /GZ -Zm200 /c
# ADD BASE RSC /l 0x409
# ADD RSC /l 0x409 /d "_DEBUG"
BSC32=bscmake.exe
# ADD BSC32 /nologo
LIB32=link.exe -lib
# ADD LIB32 /nologo

!ENDIF 

# Begin Target

# Name "contour - Win32 Release"
# Name "contour - Win32 Debug"
# Begin Group "Source Files"

# PROP Default_Filter "cpp;c;cxx;rc;def;r;odl;idl;hpj;bat"
# Begin Source File

SOURCE=bin.cpp
# End Source File
# Begin Source File

SOURCE=bucketsearch.cpp
# End Source File
# Begin Source File

SOURCE=cellqueue.cpp
# End Source File
# Begin Source File

SOURCE=conplot.cpp
# End Source File
# Begin Source File

SOURCE=conplot2d.cpp
# End Source File
# Begin Source File

SOURCE=conplot3d.cpp
# End Source File
# Begin Source File

SOURCE=conplot_p.cpp
# End Source File
# Begin Source File

SOURCE=conplotreg2.cpp
# End Source File
# Begin Source File

SOURCE=conplotreg3.cpp
# End Source File
# Begin Source File

SOURCE=contour.cpp
# End Source File
# Begin Source File

SOURCE=contour2d.cpp
# End Source File
# Begin Source File

SOURCE=contour3d.cpp
# End Source File
# Begin Source File

SOURCE=cubes.cpp
# End Source File
# Begin Source File

SOURCE=data.cpp
# End Source File
# Begin Source File

SOURCE=datareg2.cpp
# End Source File
# Begin Source File

SOURCE=datareg3.cpp
# End Source File
# Begin Source File

SOURCE=dataslc.cpp
# End Source File
# Begin Source File

SOURCE=datavol.cpp
# End Source File
# Begin Source File

SOURCE=dict.c
# End Source File
# Begin Source File

SOURCE=dirseeds.cpp
# End Source File
# Begin Source File

SOURCE=dirseedsreg2.cpp
# End Source File
# Begin Source File

SOURCE=edgehash.cpp
# End Source File
# Begin Source File

SOURCE=hash.cpp
# End Source File
# Begin Source File

SOURCE=inttree.cpp
# End Source File
# Begin Source File

SOURCE=ipqueue.cpp
# End Source File
# Begin Source File

SOURCE=iqueue.cpp
# End Source File
# Begin Source File

SOURCE=queue.cpp
# End Source File
# Begin Source File

SOURCE=range.cpp
# End Source File
# Begin Source File

SOURCE=rangeprop.cpp
# End Source File
# Begin Source File

SOURCE=rangesweep.cpp
# End Source File
# Begin Source File

SOURCE=regprop.cpp
# End Source File
# Begin Source File

SOURCE=regprop2.cpp
# End Source File
# Begin Source File

SOURCE=respprop2.cpp
# End Source File
# Begin Source File

SOURCE=seedall.cpp
# End Source File
# Begin Source File

SOURCE=seedcells.cpp
# End Source File
# Begin Source File

SOURCE=seedchkr2.cpp
# End Source File
# Begin Source File

SOURCE=seedchkr3.cpp
# End Source File
# Begin Source File

SOURCE=seeddirreg3.cpp
# End Source File
# Begin Source File

SOURCE=segtree.cpp
# End Source File
# Begin Source File

SOURCE=shelf.cpp
# End Source File
# Begin Source File

SOURCE=spectrumtest.cpp
# End Source File
# Begin Source File

SOURCE=squeue.cpp
# End Source File
# End Group
# Begin Group "Header Files"

# PROP Default_Filter "h;hpp;hxx;hm;inl"
# Begin Source File

SOURCE=basic.h
# End Source File
# Begin Source File

SOURCE=bin.h
# End Source File
# Begin Source File

SOURCE=bucketsearch.h
# End Source File
# Begin Source File

SOURCE=cellqueue.h
# End Source File
# Begin Source File

SOURCE=cellsearch.h
# End Source File
# Begin Source File

SOURCE=commdata.h
# End Source File
# Begin Source File

SOURCE=compute.h
# End Source File
# Begin Source File

SOURCE=conplot.h
# End Source File
# Begin Source File

SOURCE=conplot2d.h
# End Source File
# Begin Source File

SOURCE=conplot3d.h
# End Source File
# Begin Source File

SOURCE=conplotreg2.h
# End Source File
# Begin Source File

SOURCE=conplotreg3.h
# End Source File
# Begin Source File

SOURCE=contour.h
# End Source File
# Begin Source File

SOURCE=contour2d.h
# End Source File
# Begin Source File

SOURCE=contour3d.h
# End Source File
# Begin Source File

SOURCE=cubes.h
# End Source File
# Begin Source File

SOURCE=data.h
# End Source File
# Begin Source File

SOURCE=datareg2.h
# End Source File
# Begin Source File

SOURCE=datareg3.h
# End Source File
# Begin Source File

SOURCE=dataset.h
# End Source File
# Begin Source File

SOURCE=datasetreg2.h
# End Source File
# Begin Source File

SOURCE=datasetreg3.h
# End Source File
# Begin Source File

SOURCE=datasetslc.h
# End Source File
# Begin Source File

SOURCE=datasetvol.h
# End Source File
# Begin Source File

SOURCE=dataslc.h
# End Source File
# Begin Source File

SOURCE=datavol.h
# End Source File
# Begin Source File

SOURCE=dict.h
# End Source File
# Begin Source File

SOURCE=dirseeds.h
# End Source File
# Begin Source File

SOURCE=dirseedsreg2.h
# End Source File
# Begin Source File

SOURCE=edgehash.h
# End Source File
# Begin Source File

SOURCE=endian_io.h
# End Source File
# Begin Source File

SOURCE=hash.h
# End Source File
# Begin Source File

SOURCE=inttree.h
# End Source File
# Begin Source File

SOURCE=ipqueue.h
# End Source File
# Begin Source File

SOURCE=iqueue.h
# End Source File
# Begin Source File

SOURCE=queue.h
# End Source File
# Begin Source File

SOURCE=range.h
# End Source File
# Begin Source File

SOURCE=rangeprop.h
# End Source File
# Begin Source File

SOURCE=rangesweep.h
# End Source File
# Begin Source File

SOURCE=regprop.h
# End Source File
# Begin Source File

SOURCE=regprop2.h
# End Source File
# Begin Source File

SOURCE=respprop2.h
# End Source File
# Begin Source File

SOURCE=seedall.h
# End Source File
# Begin Source File

SOURCE=seedcells.h
# End Source File
# Begin Source File

SOURCE=seedchkr2.h
# End Source File
# Begin Source File

SOURCE=seedchkr3.h
# End Source File
# Begin Source File

SOURCE=seeddirreg3.h
# End Source File
# Begin Source File

SOURCE=segtree.h
# End Source File
# Begin Source File

SOURCE=shelf.h
# End Source File
# Begin Source File

SOURCE=squeue.h
# End Source File
# Begin Source File

SOURCE=utilities.h
# End Source File
# Begin Source File

SOURCE=vtkmarchingcubescases.h
# End Source File
# End Group
# Begin Group "Resource Files"

# PROP Default_Filter "ico;cur;bmp;dlg;rc2;rct;bin;rgs;gif;jpg;jpeg;jpe"
# End Group
# Begin Group "Generated"

# PROP Default_Filter "moc"
# End Group
# End Target
# End Project
