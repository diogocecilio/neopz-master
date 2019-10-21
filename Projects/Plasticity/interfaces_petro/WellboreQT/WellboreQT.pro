#-------------------------------------------------
#
# Project created by QtCreator 2013-11-28T10:41:54
#
#-------------------------------------------------
#Begin VTK
#INCLUDEPATH += "C://Program Files (x86)//VTK//include//vtk-5.10"
#LIBS += -L"C:/Program Files (x86)/VTK/lib/vtk-5.10" -lvtkCommon -lvtksys -lQVTK -lvtkRendering -lvtkGraphics -lvtkIO -lvtkInfovis -lvtkViews -lvtkFiltering -lvtkHybrid -lvtkWidgets
INCLUDEPATH += "/usr/include/vtk-5.8"
LIBS += -L"/usr/lib/vtk-5.8" -lvtkCommon -lvtksys -lQVTK -lvtkRendering -lvtkGraphics -lvtkIO -lvtkInfovis -lvtkViews -lvtkFiltering -lvtkHybrid -lvtkWidgets
#End VTK

DEFINES = DEBUG

#Begin PZ
INCLUDEPATH += "/dados/GOOGLE_PZ/PZLIB/include"
LIBS += -L"/dados/GOOGLE_PZ/PZLIB/" /dados/GOOGLE_PZ/PZLIB/libpz_dbg.a
#LIBS += -L"/dados/GOOGLE_PZ/PZLIB/" /dados/GOOGLE_PZ/PZLIB/libpz_rel.a
DEFINES += PZSOURCEDIR=\"/dados/GOOGLE_PZ/neopz_build_teste\" REFPATTERNDIR=\"/dados/GOOGLE_PZ/neopz_build_teste\"
DEFINES += USING_BOOST _AUTODIFF USING_FAD USING_METIS BUILD_UNITTESTING BUILD_TUTORIAL BUILD_PLASTICITY_MATERIALS REALdouble STATEdouble
#End PZ

#Begin LOG4CXX
#DEFINES += LOG4CXX
INCLUDEPATH += "/usr/local/include/log4cxx"
LIBS += -L"/usr/local/lib/" -llog4cxx
#End LOG4CXX

#Begin METIS
INCLUDEPATH += "/usr/local/include/"
LIBS += -L"/usr/local/lib/" -lmetis
#End METIS

#Begin FAD
INCLUDEPATH += "/dados/GOOGLE_PZ/externallibs/FAD/"
INCLUDEPATH += "/dados/GOOGLE_PZ/externallibs/FAD/Fad"
INCLUDEPATH += "/dados/GOOGLE_PZ/externallibs/FAD/TinyFad"
INCLUDEPATH += "/dados/GOOGLE_PZ/externallibs/FAD/TinyFadET"
#End FAD

#Begin PThread
#INCLUDEPATH += "C:/Users/Raul/Downloads/externallibs/include"
#LIBS += -L"C:/Users/Raul/Downloads/externallibs/lib" -lpthread
INCLUDEPATH += "/usr/local/include"
LIBS += -L"/usr/local/lib/" -lPThread
#End PThread



#PLASTICITY FILES
INCLUDEPATH += "/dados/GOOGLE_PZ/neopz_build_teste/Projects/Plasticity"


QT       += core gui

TARGET = WellboreQT
TEMPLATE = app

TRANSLATIONS = tr_ptbr.ts


SOURCES += main.cpp\
        mainwindow.cpp \
    evaluatedw.cpp \
    parametersdw.cpp \
    loadparamsdw.cpp \
    geometryparamsdw.cpp \
    identifyellipsdw.cpp \
/dados/GOOGLE_PZ/neopz_build_teste/Projects/Plasticity/WellBoreAnalysis.cpp \
/dados/GOOGLE_PZ/neopz_build_teste/Projects/Plasticity/TPBrBiotForce.cpp \
/dados/GOOGLE_PZ/neopz_build_teste/Projects/Plasticity/GeoMeshClass.cpp \
/dados/GOOGLE_PZ/neopz_build_teste/Projects/Plasticity/pzelastoplasticSest2D.cpp \
/dados/GOOGLE_PZ/neopz_build_teste/Projects/Plasticity/pzelasticSest2D.cpp \
    wellboreinteractorstyle.cpp \
    pwbdw.cpp \
    initialconfigdialog.cpp \
    historydw.cpp \
    prefinedw.cpp \
    hrefinedw.cpp \
    fluiddw.cpp \
    poromechdw.cpp \
    meshdw.cpp \
    tpbrunitinput.cpp \
    externallibs/units-2.11/units.c \
    externallibs/units-2.11/parse.tab.c \
    acidjobdw.cpp \
    casingdw.cpp \
    solutiondw.cpp

HEADERS  += mainwindow.h \
    evaluatedw.h \
    parametersdw.h \
    loadparamsdw.h \
    geometryparamsdw.h \
    identifyellipsdw.h \
/dados/GOOGLE_PZ/neopz_build_teste/Projects/Plasticity/WellBoreAnalysis.h \
/dados/GOOGLE_PZ/neopz_build_teste/Projects/Plasticity/TPBrBiotForce.h \
/dados/GOOGLE_PZ/neopz_build_teste/Projects/Plasticity/GeoMeshClass.h \
/dados/GOOGLE_PZ/neopz_build_teste/Projects/Plasticity/pzelastoplasticSest2D.h \
/dados/GOOGLE_PZ/neopz_build_teste/Projects/Plasticity/pzelasticSest2D.h \
    wellboreinteractorstyle.h \
    pwbdw.h \
    initialconfigdialog.h \
    historydw.h \
    prefinedw.h \
    hrefinedw.h \
    fluiddw.h \
    poromechdw.h \
    meshdw.h \
    greekletters.h \
    tpbrunitinput.h \
    externallibs/units-2.11/units.h \
    acidjobdw.h \
    casingdw.h \
    solutiondw.h

FORMS    += mainwindow.ui \
    evaluatedw.ui \
    parametersdw.ui \
    geometryparamsdw.ui \
    loadparamsdw.ui \
    identifyellipsdw.ui \
    pwbdw.ui \
    initialconfigdialog.ui \
    historydw.ui \
    prefinedw.ui \
    hrefinedw.ui \
    fluiddw.ui \
    poromechdw.ui \
    meshdw.ui \
    acidjobdw.ui \
    casingdw.ui \
    solutiondw.ui

RESOURCES += \
    SEST2D.qrc

