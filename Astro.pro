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
    CortexLikelihood.cpp \
    Splines.cpp

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
    MCMC.h \
    Splines.h

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
    run/p_013/o_013 \
    run/p_013/parameters_013.txt \
    run/p_114/parameters_114.txt \
    run/p_012/o_012 \
    run/p_012/parameters_012.txt \
    run/p_111/parameters_111.txt \
    run/p_111/o_111 \
    run/p_112/o_112 \
    run/p_114/o_114 \
    run/opt_script \
    run/p_011/o_011 \
    run/p_011/parameters_011.txt \
    run/p_000/parameters_000.txt \
    run/p_000/o_000 \
    run/experiment_0_1 \
    run/p_021/parameters_021.txt \
    run/p_022/parameters_022.txt \
    run/p_023/parameters_023.txt \
    run/p_121/parameters_121.txt \
    run/p_123/parameters_123.txt \
    run/p_122/parameters_122.txt \
    run/p_124/parameters_124.txt \
    run/p_031/parameters_031.txt \
    run/p_131/parameters_131.txt \
    run/p_132/parameters_132.txt \
    run/p_141/parameters_141.txt \
    run/3dplCL.txt \
    run/p_114_24/parameters_114024.txt \
    run/p_112_22/parameters_112_22.txt \
    run/experiment_log_250_1 \
    run/experiment_log_50_1 \
    run/experiment_log_100_1 \
    run/experiment_log_25_1 \
    run/p_144/parameters_144.txt \
    run/p_114_24_44/parameters_114_24_44.txt \
    run/p_142/parameters_142.txt \
    run/p_013_23_31/parameters_013_23_31.txt \
    run/p_012_22/parameters_012_22.txt \
    run/p_114_24_32_44/parameters_114_24_32_44.txt \
    run/p_012_22/o_012_22 \
    run/p_013_23_31/o_013_23_31 \
    run/p_021/o_021 \
    run/p_022/o_022 \
    run/p_023/o_023 \
    run/p_031/o_031 \
    run/p_112_22/o_112_22 \
    run/p_114_24/o_114_24 \
    run/p_114_24_32_44/o_114_24_32_44 \
    run/p_114_24_44/o_114_24_44 \
    run/p_121/o_121 \
    run/p_122/o_122 \
    run/p_123/o_123 \
    run/p_124/o_124 \
    run/p_131/o_131 \
    run/p_132/o_132 \
    run/p_141/o_141 \
    run/p_142/o_142 \
    run/p_144/o_144 \
    run/runModel_0.sh \
    run/runModel_1.sh \
    run/p_013_23/parameters_013_23.txt \
    run/p_013_23/o_13_23

win32{
LIBS += -L$$PWD/bin -lcygblas \
        -L$$PWD/bin -lcyglapack
} else {
LIBS += -L$$PWD/bin -lblas \
        -L$$PWD/bin -llapack
}
