TEMPLATE = app
CONFIG += console
CONFIG -= app_bundle
CONFIG -= qt

QMAKE_CXXFLAGS_RELEASE += -std=c++11 #-lpthread
QMAKE_CXXFLAGS_DEBUG += -std=c++11 #-lpthread

QMAKE_CXX = g++-5


SOURCES += main.cpp \
    CommandManager.cpp \
    Models.cpp \
    CortexSimulation.cpp \
    CortexMeasure.cpp \
    Parameters.cpp \
    MatrixInverse.cpp \
    matrixCholesky.cpp \
    LevenbergMarquardt.cpp \
    BayesIteration.cpp \
    CortexLikelihood.cpp

#include(deployment.pri)
#qtcAddDeployment()

HEADERS += \
    CommandManager.h \
    Models.h \
    CortexSimulation.h \
    CortexMeasure.h \
    Parameters.h \
    MatrixInverse.h \
    LevenbergMarquardt.h \
    BayesIteration.h \
    CortexLikelihood.h \
    BaseClass.h \
    MCMC.h

DISTFILES += \
    run/script \
    run/3dpl.txt \
    run/minimal_prior.txt \
    run/7dpl.txt \
    run/sim_zero.txt \
    run/simulate_script \
    run/parameters_zero.txt \
    run/write_script \
    run/parameters_10 \
    run/opt_script \
    run/parameters_10.txt \
    run/opt_modelzero \
    run/opt_model10 \
    run/parameters_011.txt \
    run/parameters_012.txt \
    run/parameters_013.txt \
    run/.txt \
    run/parameters_013.txt \
    run/parameters_111.txt \
    run/parameters_112.txt \
    run/parameters_113.txt \
    run/parameters_114.txt \
    run/parameters_116.txt \
    Makefile \
    run/7dpl2.txt \
    run/o.sh

win32{
LIBS += -L$$PWD/bin -lcygblas \
        -L$$PWD/bin -lcyglapack
} else {
LIBS += -L$$PWD/bin -lblas \
        -L$$PWD/bin -llapack
}
