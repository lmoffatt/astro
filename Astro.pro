TEMPLATE = app
CONFIG += console
CONFIG -= app_bundle
CONFIG -= qt

QMAKE_CXXFLAGS_RELEASE += -std=c++11 #-lpthread
QMAKE_CXXFLAGS_DEBUG += -std=c++11 #-lpthread


SOURCES += main.cpp \
    CortexState.cpp \
    Script.cpp \
    CommandManager.cpp \
    read.cpp \
    AlignCommand.cpp

include(deployment.pri)
qtcAddDeployment()

HEADERS += \
    Astrocyte.h \
    CortexState.h \
    Script.h \
    CommandManager.h \
    read.h \
    AlignCommand.h

DISTFILES += \
    run/script \
    run/3dpl.txt \
    run/minimal_prior.txt \
    run/7dpl.txt \
    run/sim_zero.txt \
    run/simulate_script \
    run/parameters_zero.txt

