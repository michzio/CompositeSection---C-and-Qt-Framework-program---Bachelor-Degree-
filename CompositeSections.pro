#-------------------------------------------------
#
# Project created by QtCreator 2013-07-27T09:41:13
#
#-------------------------------------------------

QT       += core gui
CONFIG   += c++11
QT       += opengl

greaterThan(QT_MAJOR_VERSION, 4): QT += widgets
                                  QT += xml

TARGET = CompositeSections
TEMPLATE = app


SOURCES += main.cpp\
        mainwindow.cpp \
    points.cpp \
    material.cpp \
    surface.cpp \
    line.cpp \
    fibergroup.cpp \
    point.cpp \
    section.cpp \
    sectionsolver.cpp \
    vector.cpp \
    arc.cpp \
    concrete.cpp \
    reinforcement.cpp \
    structuralsteel.cpp \
    polygonintersection.cpp \
    gpc.c \
    combinationgenerator.cpp \
    mouseeventshandler.cpp \
    mygraphicsview.cpp \
    strainprofile.cpp \
    segment.cpp \
    gausslegendre.cpp \
    solutionwindow.cpp \
    glwidget.cpp \
    vertex3d.cpp \
    gridgraphicsscene.cpp \
    materialwindow.cpp \
    fafitisconcrete.cpp \
    fafitisreinforcement.cpp \
    hollow.cpp \
    frp.cpp

HEADERS  += mainwindow.h \
    ../material.h \
    points.h \
    material.h \
    surface.h \
    line.h \
    fibergroup.h \
    point.h \
    section.h \
    sectionsolver.h \
    vector.h \
    arc.h \
    concrete.h \
    reinforcement.h \
    structuralsteel.h \
    polygonintersection.h \
    gpc.h \
    combinationgenerator.h \
    mouseeventshandler.h \
    mygraphicsview.h \
    strainprofile.h \
    segment.h \
    gausslegendre.h \
    solutionwindow.h \
    glwidget.h \
    vertex3d.h \
    flags.h \
    gridgraphicsscene.h \
    materialwindow.h \
    fafitisconcrete.h \
    fafitisreinforcement.h \
    hollow.h \
    frp.h

FORMS    += mainwindow.ui \
    solutionwindow.ui \
    materialwindow.ui
