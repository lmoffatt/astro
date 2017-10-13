TEMPLATE = app
CONFIG += console
CONFIG -= app_bundle
CONFIG -= qt

CONFIG += c++14

# remove possible other optimization flags
QMAKE_CXXFLAGS_RELEASE -= -O
QMAKE_CXXFLAGS_RELEASE -= -O1
QMAKE_CXXFLAGS_RELEASE -= -O2

# add the desired -O3 if not present
QMAKE_CXXFLAGS_RELEASE += -O2 -Werror
QMAKE_CXXFLAGS_RELEASE += -lpthread
QMAKE_CXXFLAGS_DEBUG += -lpthread


QMAKE_CXXFLAGS += -fopenmp
QMAKE_LFLAGS +=  -fopenmp
LIBS += -lpthread

QMAKE_CXXFLAGS_RELEASE+="-D__STRICT_ANSI__"
QMAKE_CXXFLAGS_DEBUG+="-D__STRICT_ANSI__"

#QMAKE_CXX = g++


SOURCES += main.cpp \
    CommandManager.cpp \
    Models.cpp \
    CortexSimulation.cpp \
    CortexMeasure.cpp \
    Parameters.cpp \
    MatrixInverse.cpp \
    matrixCholesky.cpp \
    LevenbergMarquardt.cpp \
    CortexLikelihood.cpp \
    Splines.cpp \
    Evidence.cpp \
    Matrix.cpp \
    Distributions.cpp

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
    CortexLikelihood.h \
    BaseClass.h \
    Splines.h \
    Evidence.h \
    Matrix.h \
    Distributions.h \
    Optimization.h \
    Optimization_BFGS.h \
    mySerializer.h \
    myOutputSerializer.h \
    myInputSerializer.h \
    myOrderOperators.h \
    myTuples.h \
    myCommandManagement.h \
    RungeKutta.h \
    MatrixBanded.h \
    MatrixBlock.h

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
    run/p_113/parameters_113.txt \
    run/p_115/parameters_115.txt \
    run/p_113/o_113 \
    run/p_115/o_115 \
    run/runModelNew.sh \
    run/p_114_24/parameters_114_24.txt \
    run/p_115_22/parameters_115_22.txt \
    run/p_115_22/o_115_22 \
    run/p_125/parameters_125.txt \
    run/p_125/o_125 \
    run/p_115_25/parameters_115_25.txt \
    run/p_115_25/o_115_25 \
    run/p_113_42/parameters_113_42 \
    run/p_112_22_31/parameters_112_22_31 \
    run/p_112_22_31/o_112_22_31 \
    run/p_113_42/o_113_42 \
    run/p_013_23/o_013_23 \
    run/p_051/parameters_051.txt \
    run/p_012_51/parameters_012_51.txt \
    run/p_013_51/parameters_013_51.txt \
    run/p_051/o_051 \
    run/p_013_51/o_013_51 \
    run/p_012_51/o_012_51 \
    run/p_151/parameters_151.txt \
    run/p_112_51/parameters_112_51.txt \
    run/p_200/parameters_200.txt \
    run/p_212/parameters_212.txt \
    run/p_214/parameters_214.txt \
    run/p_200/o_200 \
    run/p_212/o_212 \
    run/p_214/o_214 \
    run/p_151/o_151 \
    run/P_152/o_152 \
    run/p_112_51/o_112_51 \
    run/p_112_52/parameters_112_52.txt \
    run/p_112_52/o_112_52 \
    run/p_213/parameters_213.txt \
    run/p_215/parameters_215.txt \
    run/p_211/parameters_211.txt \
    run/p_211/o_211 \
    run/p_213/o_213 \
    run/p_215/o_215 \
    run/runModel_2.sh \
    run/p_152/parameters_152.txt \
    run/p_152/o_152 \
    run/p_213_52/parameters_213_52.txt \
    run/p_213_51/parameters_213_51 \
    run/p_213_23/parameters_213_23 \
    run/local_script \
    run/mcmc_script \
    debug/Makefile \
    run/7dpl_resultados.txt \
    run/p_000m/parameters_000m.txt \
    run/p_100m/parameters_100m.txt \
    run/New_experiment_log_250_1 \
    run/evidence_script \
    run/runEvidence_0.sh \
    run/DATOS_7dpl_resultados.txt \
    run/DATOS_3dpl_CONTROL_resultados.txt \
    run/DATOS_3dpl_resultados.txt \
    run/E_250_50 \
    run/E_250 \
    run/p_000m/o_000m \
    run/p_100m/o_100m \
    run/p_114_24_32_44m/o_114_24_32_44m \
    run/p_114_24_32_44m/parameters_114_24_32_44m.txt \
    run/p_013_23_31m/o_013_23_31m \
    run/p_013_23_31m/parameters_013_23_31m.txt \
    run/m10/model_10.txt \
    run/m13/model_13.txt \
    run/m01/model_01.txt \
    run/m12/model_12.txt \
    run/m11/model_11.txt \
    run/m101/model_101.txt \
    run/m03/model_03.txt \
    run/m02/model_02.txt \
    run/m01/o_m01 \
    run/m02/o_m02 \
    run/m03/o_m03 \
    run/m10/o_m10 \
    run/m11/o_m11 \
    run/m13/o_m13 \
    run/m12/o_m12 \
    run/evidence_script_2 \
    run/evidence_script_3 \
    run/evidence_script_4 \
    run/0 \
    run/evidence_script_1 \
    run/evidence_script_5 \
    run/evidence_script_6 \
    run/evidence_script_7 \
    run/evidence_script \
    run/m01/o_m01_1 \
    run/m01/o_m01_2 \
    run/m01/o_m01_3 \
    run/0 \
    run/runLandas_0.sh \
    run/m01/o_m01_4 \
    run/m01/o_m01_5 \
    run/m01/o_m01_6 \
    run/m01/o_m01_7 \
    run/paralel/Makefile \
    run/m01/p_m01 \
    run/m01/p_m02 \
    run/m01/p_m03 \
    run/m10/p_m10 \
    run/m101/p_m101 \
    run/runParallel_0.sh \
    run/m01/p_m01_1 \
    run/m10/p_m10_1 \
    run/m101/p_m101_1 \
    run/m01/p_m01_2 \
    run/m10/p_m10_2 \
    run/m101/p_m101_2 \
    run/m03/p_m03 \
    run/m13/p_m13 \
    run/m02/p_m02 \
    run/m11/p_m11 \
    run/m12/p_m12 \
    run/E_250_50_1 \
    run/m02/p_m02 \
    run/m01/s_m01 \
    slurm/Makefile \
    run/m02/s_m02

win32{
LIBS += -L$$PWD/bin -lcygblas \
        -L$$PWD/bin -lcyglapack
} else {
LIBS += -L$$PWD/bin -lblas \
        -L$$PWD/bin -llapack
}
