TEMPLATE = app
CONFIG += console
CONFIG -= app_bundle
CONFIG -= qt

QMAKE_CXXFLAGS_RELEASE += -std=c++11 #-lpthread
QMAKE_CXXFLAGS_DEBUG += -std=c++11 #-lpthread

QMAKE_CXX = g++


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
    Makefile \
    run/7dpl2.txt \
    run/p_100/o_100 \
    run/p_100/parameters_10.txt \
    run/p_011/parameters_111.txt \
    run/p_112/parameters_112.txt \
    run/p_116/parameters_116.txt \
    run/p_116/o_116 \
    run/p_013/o_013 \
    run/p_013/parameters_013.txt \
    run/p_114/parameters_114.txt \
    run/p_012/o_012 \
    run/p_012/parameters_012.txt \
    run/p_111/parameters_111.txt \
    run/p_111/o_111 \
    run/p_112/o_112 \
    run/p_114/o_114 \
    run/p_000/o_011 \
    run/p_000/parameters_011.txt \
    run/p_011/parameters_000.txt \
    run/p_011/o_000

win32{
LIBS += -L$$PWD/bin -lcygblas \
        -L$$PWD/bin -lcyglapack
} else {
LIBS += -L$$PWD/bin -lblas \
        -L$$PWD/bin -llapack
}
