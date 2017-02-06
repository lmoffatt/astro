#############################################################################
# Makefile for building: Astro
# Generated by qmake (3.0) (Qt 5.7.0)
# Project:  ../Astro/Astro.pro
# Template: app
# Command: /home/luciano/Qt5.7.0/5.7/gcc_64/bin/qmake -spec linux-g++ -o Makefile ../Astro/Astro.pro
#############################################################################

MAKEFILE      = Makefile

####### Compiler, tools and options

CC            = gcc
CXX           = g++
DEFINES       =
CFLAGS        = -pipe -O2 -Wall -W -fPIC $(DEFINES)
CXXFLAGS      = -pipe -fopenmp -O2 -std=c++11 -std=gnu++11 -Wall -W -fPIC $(DEFINES)
INCPATH       = -I../Astro -I. -I../../../../Qt5.7.0/5.7/gcc_64/mkspecs/linux-g++
QMAKE         = /home/luciano/Qt5.7.0/5.7/gcc_64/bin/qmake
DEL_FILE      = rm -f
CHK_DIR_EXISTS= test -d
MKDIR         = mkdir -p
COPY          = cp -f
COPY_FILE     = cp -f
COPY_DIR      = cp -f -R
INSTALL_FILE  = install -m 644 -p
INSTALL_PROGRAM = install -m 755 -p
INSTALL_DIR   = cp -f -R
DEL_FILE      = rm -f
SYMLINK       = ln -f -s
DEL_DIR       = rmdir
MOVE          = mv -f
TAR           = tar -cf
COMPRESS      = gzip -9f
DISTNAME      = Astro1.0.0
DISTDIR = /home/luciano/Data/celulas/Astro/build-Astro-gcc_6-Release/.tmp/Astro1.0.0
LINK          = g++
LFLAGS        = -fopenmp -Wl,-O1
LIBS          = $(SUBLIBS) -lgomp -lpthread -L/home/luciano/Data/celulas/Astro/Astro/bin -lblas -llapack
AR            = ar cqs
RANLIB        =
SED           = sed
STRIP         = strip

####### Output directory

OBJECTS_DIR   = ./

####### Files

SOURCES       = ../Astro/main.cpp \
../Astro/CommandManager.cpp \
../Astro/Models.cpp \
		../Astro/CortexSimulation.cpp \
		../Astro/CortexMeasure.cpp \
		../Astro/Parameters.cpp \
		../Astro/MatrixInverse.cpp \
		../Astro/matrixCholesky.cpp \
		../Astro/LevenbergMarquardt.cpp \
		../Astro/CortexLikelihood.cpp \
		../Astro/Splines.cpp \
		../Astro/Evidence.cpp \
		../Astro/Matrix.cpp \
		../Astro/Distributions.cpp
OBJECTS       = main.o \
CommandManager.o \
Models.o \
		CortexSimulation.o \
		CortexMeasure.o \
		Parameters.o \
		MatrixInverse.o \
		matrixCholesky.o \
		LevenbergMarquardt.o \
		CortexLikelihood.o \
		Splines.o \
		Evidence.o \
		Matrix.o \
		Distributions.o
DIST          = ../Astro/run/script \
../Astro/run/3dpl.txt \
../Astro/run/minimal_prior.txt \
		../Astro/run/7dpl.txt \
		../Astro/run/sim_zero.txt \
		../Astro/run/simulate_script \
		../Astro/run/parameters_zero.txt \
		../Astro/run/write_script \
		../Astro/run/parameters_10 \
		../Astro/run/opt_script \
		../Astro/run/parameters_10.txt \
		../Astro/run/opt_modelzero \
		../Astro/run/opt_model10 \
		../Astro/Makefile \
		../Astro/run/7dpl2.txt \
		../Astro/run/p_100/o_100 \
		../Astro/run/p_100/parameters_10.txt \
		../Astro/run/p_011/parameters_111.txt \
		../Astro/run/p_112/parameters_112.txt \
		../Astro/run/p_013/o_013 \
		../Astro/run/p_013/parameters_013.txt \
		../Astro/run/p_114/parameters_114.txt \
		../Astro/run/p_012/o_012 \
		../Astro/run/p_012/parameters_012.txt \
		../Astro/run/p_111/parameters_111.txt \
		../Astro/run/p_111/o_111 \
		../Astro/run/p_112/o_112 \
		../Astro/run/p_114/o_114 \
		../Astro/run/opt_script \
		../Astro/run/p_011/o_011 \
		../Astro/run/p_011/parameters_011.txt \
		../Astro/run/p_000/parameters_000.txt \
		../Astro/run/p_000/o_000 \
		../Astro/run/experiment_0_1 \
		../Astro/run/p_021/parameters_021.txt \
		../Astro/run/p_022/parameters_022.txt \
		../Astro/run/p_023/parameters_023.txt \
		../Astro/run/p_121/parameters_121.txt \
		../Astro/run/p_123/parameters_123.txt \
		../Astro/run/p_122/parameters_122.txt \
		../Astro/run/p_124/parameters_124.txt \
		../Astro/run/p_031/parameters_031.txt \
		../Astro/run/p_131/parameters_131.txt \
		../Astro/run/p_132/parameters_132.txt \
		../Astro/run/p_141/parameters_141.txt \
		../Astro/run/3dplCL.txt \
		../Astro/run/p_112_22/parameters_112_22.txt \
		../Astro/run/experiment_log_250_1 \
		../Astro/run/experiment_log_50_1 \
		../Astro/run/experiment_log_100_1 \
		../Astro/run/experiment_log_25_1 \
		../Astro/run/p_144/parameters_144.txt \
		../Astro/run/p_114_24_44/parameters_114_24_44.txt \
		../Astro/run/p_142/parameters_142.txt \
		../Astro/run/p_013_23_31/parameters_013_23_31.txt \
		../Astro/run/p_012_22/parameters_012_22.txt \
		../Astro/run/p_114_24_32_44/parameters_114_24_32_44.txt \
		../Astro/run/p_012_22/o_012_22 \
		../Astro/run/p_013_23_31/o_013_23_31 \
		../Astro/run/p_021/o_021 \
		../Astro/run/p_022/o_022 \
		../Astro/run/p_023/o_023 \
		../Astro/run/p_031/o_031 \
		../Astro/run/p_112_22/o_112_22 \
		../Astro/run/p_114_24/o_114_24 \
		../Astro/run/p_114_24_32_44/o_114_24_32_44 \
		../Astro/run/p_114_24_44/o_114_24_44 \
		../Astro/run/p_121/o_121 \
		../Astro/run/p_122/o_122 \
		../Astro/run/p_123/o_123 \
		../Astro/run/p_124/o_124 \
		../Astro/run/p_131/o_131 \
		../Astro/run/p_132/o_132 \
		../Astro/run/p_141/o_141 \
		../Astro/run/p_142/o_142 \
		../Astro/run/p_144/o_144 \
		../Astro/run/runModel_0.sh \
		../Astro/run/runModel_1.sh \
		../Astro/run/p_013_23/parameters_013_23.txt \
		../Astro/run/p_113/parameters_113.txt \
		../Astro/run/p_115/parameters_115.txt \
		../Astro/run/p_113/o_113 \
		../Astro/run/p_115/o_115 \
		../Astro/run/runModelNew.sh \
		../Astro/run/p_114_24/parameters_114_24.txt \
		../Astro/run/p_115_22/parameters_115_22.txt \
		../Astro/run/p_115_22/o_115_22 \
		../Astro/run/p_125/parameters_125.txt \
		../Astro/run/p_125/o_125 \
		../Astro/run/p_115_25/parameters_115_25.txt \
		../Astro/run/p_115_25/o_115_25 \
		../Astro/run/p_113_42/parameters_113_42 \
		../Astro/run/p_112_22_31/parameters_112_22_31 \
		../Astro/run/p_112_22_31/o_112_22_31 \
		../Astro/run/p_113_42/o_113_42 \
		../Astro/run/p_013_23/o_013_23 \
		../Astro/run/p_051/parameters_051.txt \
		../Astro/run/p_012_51/parameters_012_51.txt \
		../Astro/run/p_013_51/parameters_013_51.txt \
		../Astro/run/p_051/o_051 \
		../Astro/run/p_013_51/o_013_51 \
		../Astro/run/p_012_51/o_012_51 \
		../Astro/run/p_151/parameters_151.txt \
		../Astro/run/p_112_51/parameters_112_51.txt \
		../Astro/run/p_200/parameters_200.txt \
		../Astro/run/p_212/parameters_212.txt \
		../Astro/run/p_214/parameters_214.txt \
		../Astro/run/p_200/o_200 \
		../Astro/run/p_212/o_212 \
		../Astro/run/p_214/o_214 \
		../Astro/run/p_151/o_151 \
		../Astro/run/P_152/o_152 \
		../Astro/run/p_112_51/o_112_51 \
		../Astro/run/p_112_52/parameters_112_52.txt \
		../Astro/run/p_112_52/o_112_52 \
		../Astro/run/p_213/parameters_213.txt \
		../Astro/run/p_215/parameters_215.txt \
		../Astro/run/p_211/parameters_211.txt \
		../Astro/run/p_211/o_211 \
		../Astro/run/p_213/o_213 \
		../Astro/run/p_215/o_215 \
		../Astro/run/runModel_2.sh \
		../Astro/run/p_152/parameters_152.txt \
		../Astro/run/p_152/o_152 \
		../Astro/run/p_213_52/parameters_213_52.txt \
		../Astro/run/p_213_51/parameters_213_51 \
		../Astro/run/p_213_23/parameters_213_23 \
		../Astro/run/local_script \
		../Astro/run/mcmc_script \
		../Astro/debug/Makefile \
		../Astro/run/7dpl_resultados.txt \
		../Astro/run/p_000m/parameters_000m.txt \
		../Astro/run/p_100m/parameters_100m.txt \
		../Astro/run/p_114_24_32_44m/parameter_114_24_32_44 \
		../Astro/run/New_experiment_log_250_1 \
		../Astro/run/evidence_script \
		../../../../Qt5.7.0/5.7/gcc_64/mkspecs/features/spec_pre.prf \
		../../../../Qt5.7.0/5.7/gcc_64/mkspecs/common/unix.conf \
		../../../../Qt5.7.0/5.7/gcc_64/mkspecs/common/linux.conf \
		../../../../Qt5.7.0/5.7/gcc_64/mkspecs/common/sanitize.conf \
		../../../../Qt5.7.0/5.7/gcc_64/mkspecs/common/gcc-base.conf \
		../../../../Qt5.7.0/5.7/gcc_64/mkspecs/common/gcc-base-unix.conf \
		../../../../Qt5.7.0/5.7/gcc_64/mkspecs/common/g++-base.conf \
		../../../../Qt5.7.0/5.7/gcc_64/mkspecs/common/g++-unix.conf \
		../../../../Qt5.7.0/5.7/gcc_64/mkspecs/qconfig.pri \
		../../../../Qt5.7.0/5.7/gcc_64/mkspecs/modules/qt_lib_3dcore.pri \
		../../../../Qt5.7.0/5.7/gcc_64/mkspecs/modules/qt_lib_3dcore_private.pri \
		../../../../Qt5.7.0/5.7/gcc_64/mkspecs/modules/qt_lib_3dextras.pri \
		../../../../Qt5.7.0/5.7/gcc_64/mkspecs/modules/qt_lib_3dextras_private.pri \
		../../../../Qt5.7.0/5.7/gcc_64/mkspecs/modules/qt_lib_3dinput.pri \
		../../../../Qt5.7.0/5.7/gcc_64/mkspecs/modules/qt_lib_3dinput_private.pri \
		../../../../Qt5.7.0/5.7/gcc_64/mkspecs/modules/qt_lib_3dlogic.pri \
		../../../../Qt5.7.0/5.7/gcc_64/mkspecs/modules/qt_lib_3dlogic_private.pri \
		../../../../Qt5.7.0/5.7/gcc_64/mkspecs/modules/qt_lib_3dquick.pri \
		../../../../Qt5.7.0/5.7/gcc_64/mkspecs/modules/qt_lib_3dquick_private.pri \
		../../../../Qt5.7.0/5.7/gcc_64/mkspecs/modules/qt_lib_3dquickextras.pri \
		../../../../Qt5.7.0/5.7/gcc_64/mkspecs/modules/qt_lib_3dquickextras_private.pri \
		../../../../Qt5.7.0/5.7/gcc_64/mkspecs/modules/qt_lib_3dquickinput.pri \
		../../../../Qt5.7.0/5.7/gcc_64/mkspecs/modules/qt_lib_3dquickinput_private.pri \
		../../../../Qt5.7.0/5.7/gcc_64/mkspecs/modules/qt_lib_3dquickrender.pri \
		../../../../Qt5.7.0/5.7/gcc_64/mkspecs/modules/qt_lib_3dquickrender_private.pri \
		../../../../Qt5.7.0/5.7/gcc_64/mkspecs/modules/qt_lib_3drender.pri \
		../../../../Qt5.7.0/5.7/gcc_64/mkspecs/modules/qt_lib_3drender_private.pri \
		../../../../Qt5.7.0/5.7/gcc_64/mkspecs/modules/qt_lib_bluetooth.pri \
		../../../../Qt5.7.0/5.7/gcc_64/mkspecs/modules/qt_lib_bluetooth_private.pri \
		../../../../Qt5.7.0/5.7/gcc_64/mkspecs/modules/qt_lib_bootstrap_private.pri \
		../../../../Qt5.7.0/5.7/gcc_64/mkspecs/modules/qt_lib_charts.pri \
		../../../../Qt5.7.0/5.7/gcc_64/mkspecs/modules/qt_lib_charts_private.pri \
		../../../../Qt5.7.0/5.7/gcc_64/mkspecs/modules/qt_lib_clucene_private.pri \
		../../../../Qt5.7.0/5.7/gcc_64/mkspecs/modules/qt_lib_concurrent.pri \
		../../../../Qt5.7.0/5.7/gcc_64/mkspecs/modules/qt_lib_concurrent_private.pri \
		../../../../Qt5.7.0/5.7/gcc_64/mkspecs/modules/qt_lib_core.pri \
		../../../../Qt5.7.0/5.7/gcc_64/mkspecs/modules/qt_lib_core_private.pri \
		../../../../Qt5.7.0/5.7/gcc_64/mkspecs/modules/qt_lib_datavisualization.pri \
		../../../../Qt5.7.0/5.7/gcc_64/mkspecs/modules/qt_lib_datavisualization_private.pri \
		../../../../Qt5.7.0/5.7/gcc_64/mkspecs/modules/qt_lib_dbus.pri \
		../../../../Qt5.7.0/5.7/gcc_64/mkspecs/modules/qt_lib_dbus_private.pri \
		../../../../Qt5.7.0/5.7/gcc_64/mkspecs/modules/qt_lib_designer.pri \
		../../../../Qt5.7.0/5.7/gcc_64/mkspecs/modules/qt_lib_designer_private.pri \
		../../../../Qt5.7.0/5.7/gcc_64/mkspecs/modules/qt_lib_designercomponents_private.pri \
		../../../../Qt5.7.0/5.7/gcc_64/mkspecs/modules/qt_lib_eglfs_device_lib_private.pri \
		../../../../Qt5.7.0/5.7/gcc_64/mkspecs/modules/qt_lib_gamepad.pri \
		../../../../Qt5.7.0/5.7/gcc_64/mkspecs/modules/qt_lib_gamepad_private.pri \
		../../../../Qt5.7.0/5.7/gcc_64/mkspecs/modules/qt_lib_gui.pri \
		../../../../Qt5.7.0/5.7/gcc_64/mkspecs/modules/qt_lib_gui_private.pri \
		../../../../Qt5.7.0/5.7/gcc_64/mkspecs/modules/qt_lib_help.pri \
		../../../../Qt5.7.0/5.7/gcc_64/mkspecs/modules/qt_lib_help_private.pri \
		../../../../Qt5.7.0/5.7/gcc_64/mkspecs/modules/qt_lib_location.pri \
		../../../../Qt5.7.0/5.7/gcc_64/mkspecs/modules/qt_lib_location_private.pri \
		../../../../Qt5.7.0/5.7/gcc_64/mkspecs/modules/qt_lib_multimedia.pri \
		../../../../Qt5.7.0/5.7/gcc_64/mkspecs/modules/qt_lib_multimedia_private.pri \
		../../../../Qt5.7.0/5.7/gcc_64/mkspecs/modules/qt_lib_multimediawidgets.pri \
		../../../../Qt5.7.0/5.7/gcc_64/mkspecs/modules/qt_lib_multimediawidgets_private.pri \
		../../../../Qt5.7.0/5.7/gcc_64/mkspecs/modules/qt_lib_network.pri \
		../../../../Qt5.7.0/5.7/gcc_64/mkspecs/modules/qt_lib_network_private.pri \
		../../../../Qt5.7.0/5.7/gcc_64/mkspecs/modules/qt_lib_nfc.pri \
		../../../../Qt5.7.0/5.7/gcc_64/mkspecs/modules/qt_lib_nfc_private.pri \
		../../../../Qt5.7.0/5.7/gcc_64/mkspecs/modules/qt_lib_opengl.pri \
		../../../../Qt5.7.0/5.7/gcc_64/mkspecs/modules/qt_lib_opengl_private.pri \
		../../../../Qt5.7.0/5.7/gcc_64/mkspecs/modules/qt_lib_openglextensions.pri \
		../../../../Qt5.7.0/5.7/gcc_64/mkspecs/modules/qt_lib_openglextensions_private.pri \
		../../../../Qt5.7.0/5.7/gcc_64/mkspecs/modules/qt_lib_packetprotocol_private.pri \
		../../../../Qt5.7.0/5.7/gcc_64/mkspecs/modules/qt_lib_platformsupport_private.pri \
		../../../../Qt5.7.0/5.7/gcc_64/mkspecs/modules/qt_lib_positioning.pri \
		../../../../Qt5.7.0/5.7/gcc_64/mkspecs/modules/qt_lib_positioning_private.pri \
		../../../../Qt5.7.0/5.7/gcc_64/mkspecs/modules/qt_lib_printsupport.pri \
		../../../../Qt5.7.0/5.7/gcc_64/mkspecs/modules/qt_lib_printsupport_private.pri \
		../../../../Qt5.7.0/5.7/gcc_64/mkspecs/modules/qt_lib_purchasing.pri \
		../../../../Qt5.7.0/5.7/gcc_64/mkspecs/modules/qt_lib_purchasing_private.pri \
		../../../../Qt5.7.0/5.7/gcc_64/mkspecs/modules/qt_lib_qml.pri \
		../../../../Qt5.7.0/5.7/gcc_64/mkspecs/modules/qt_lib_qml_private.pri \
		../../../../Qt5.7.0/5.7/gcc_64/mkspecs/modules/qt_lib_qmldebug_private.pri \
		../../../../Qt5.7.0/5.7/gcc_64/mkspecs/modules/qt_lib_qmldevtools_private.pri \
		../../../../Qt5.7.0/5.7/gcc_64/mkspecs/modules/qt_lib_qmltest.pri \
		../../../../Qt5.7.0/5.7/gcc_64/mkspecs/modules/qt_lib_qmltest_private.pri \
		../../../../Qt5.7.0/5.7/gcc_64/mkspecs/modules/qt_lib_qtmultimediaquicktools_private.pri \
		../../../../Qt5.7.0/5.7/gcc_64/mkspecs/modules/qt_lib_quick.pri \
		../../../../Qt5.7.0/5.7/gcc_64/mkspecs/modules/qt_lib_quick_private.pri \
		../../../../Qt5.7.0/5.7/gcc_64/mkspecs/modules/qt_lib_quickcontrols2.pri \
		../../../../Qt5.7.0/5.7/gcc_64/mkspecs/modules/qt_lib_quickcontrols2_private.pri \
		../../../../Qt5.7.0/5.7/gcc_64/mkspecs/modules/qt_lib_quickparticles_private.pri \
		../../../../Qt5.7.0/5.7/gcc_64/mkspecs/modules/qt_lib_quicktemplates2_private.pri \
		../../../../Qt5.7.0/5.7/gcc_64/mkspecs/modules/qt_lib_quickwidgets.pri \
		../../../../Qt5.7.0/5.7/gcc_64/mkspecs/modules/qt_lib_quickwidgets_private.pri \
		../../../../Qt5.7.0/5.7/gcc_64/mkspecs/modules/qt_lib_script.pri \
		../../../../Qt5.7.0/5.7/gcc_64/mkspecs/modules/qt_lib_script_private.pri \
		../../../../Qt5.7.0/5.7/gcc_64/mkspecs/modules/qt_lib_scripttools.pri \
		../../../../Qt5.7.0/5.7/gcc_64/mkspecs/modules/qt_lib_scripttools_private.pri \
		../../../../Qt5.7.0/5.7/gcc_64/mkspecs/modules/qt_lib_scxml.pri \
		../../../../Qt5.7.0/5.7/gcc_64/mkspecs/modules/qt_lib_scxml_private.pri \
		../../../../Qt5.7.0/5.7/gcc_64/mkspecs/modules/qt_lib_sensors.pri \
		../../../../Qt5.7.0/5.7/gcc_64/mkspecs/modules/qt_lib_sensors_private.pri \
		../../../../Qt5.7.0/5.7/gcc_64/mkspecs/modules/qt_lib_serialbus.pri \
		../../../../Qt5.7.0/5.7/gcc_64/mkspecs/modules/qt_lib_serialbus_private.pri \
		../../../../Qt5.7.0/5.7/gcc_64/mkspecs/modules/qt_lib_serialport.pri \
		../../../../Qt5.7.0/5.7/gcc_64/mkspecs/modules/qt_lib_serialport_private.pri \
		../../../../Qt5.7.0/5.7/gcc_64/mkspecs/modules/qt_lib_sql.pri \
		../../../../Qt5.7.0/5.7/gcc_64/mkspecs/modules/qt_lib_sql_private.pri \
		../../../../Qt5.7.0/5.7/gcc_64/mkspecs/modules/qt_lib_svg.pri \
		../../../../Qt5.7.0/5.7/gcc_64/mkspecs/modules/qt_lib_svg_private.pri \
		../../../../Qt5.7.0/5.7/gcc_64/mkspecs/modules/qt_lib_testlib.pri \
		../../../../Qt5.7.0/5.7/gcc_64/mkspecs/modules/qt_lib_testlib_private.pri \
		../../../../Qt5.7.0/5.7/gcc_64/mkspecs/modules/qt_lib_uiplugin.pri \
		../../../../Qt5.7.0/5.7/gcc_64/mkspecs/modules/qt_lib_uitools.pri \
		../../../../Qt5.7.0/5.7/gcc_64/mkspecs/modules/qt_lib_uitools_private.pri \
		../../../../Qt5.7.0/5.7/gcc_64/mkspecs/modules/qt_lib_webchannel.pri \
		../../../../Qt5.7.0/5.7/gcc_64/mkspecs/modules/qt_lib_webchannel_private.pri \
		../../../../Qt5.7.0/5.7/gcc_64/mkspecs/modules/qt_lib_webengine.pri \
		../../../../Qt5.7.0/5.7/gcc_64/mkspecs/modules/qt_lib_webengine_private.pri \
		../../../../Qt5.7.0/5.7/gcc_64/mkspecs/modules/qt_lib_webenginecore.pri \
		../../../../Qt5.7.0/5.7/gcc_64/mkspecs/modules/qt_lib_webenginecore_private.pri \
		../../../../Qt5.7.0/5.7/gcc_64/mkspecs/modules/qt_lib_webenginecoreheaders_private.pri \
		../../../../Qt5.7.0/5.7/gcc_64/mkspecs/modules/qt_lib_webenginewidgets.pri \
		../../../../Qt5.7.0/5.7/gcc_64/mkspecs/modules/qt_lib_webenginewidgets_private.pri \
		../../../../Qt5.7.0/5.7/gcc_64/mkspecs/modules/qt_lib_websockets.pri \
		../../../../Qt5.7.0/5.7/gcc_64/mkspecs/modules/qt_lib_websockets_private.pri \
		../../../../Qt5.7.0/5.7/gcc_64/mkspecs/modules/qt_lib_webview.pri \
		../../../../Qt5.7.0/5.7/gcc_64/mkspecs/modules/qt_lib_webview_private.pri \
		../../../../Qt5.7.0/5.7/gcc_64/mkspecs/modules/qt_lib_widgets.pri \
		../../../../Qt5.7.0/5.7/gcc_64/mkspecs/modules/qt_lib_widgets_private.pri \
		../../../../Qt5.7.0/5.7/gcc_64/mkspecs/modules/qt_lib_x11extras.pri \
		../../../../Qt5.7.0/5.7/gcc_64/mkspecs/modules/qt_lib_x11extras_private.pri \
		../../../../Qt5.7.0/5.7/gcc_64/mkspecs/modules/qt_lib_xcb_qpa_lib_private.pri \
		../../../../Qt5.7.0/5.7/gcc_64/mkspecs/modules/qt_lib_xml.pri \
		../../../../Qt5.7.0/5.7/gcc_64/mkspecs/modules/qt_lib_xml_private.pri \
		../../../../Qt5.7.0/5.7/gcc_64/mkspecs/modules/qt_lib_xmlpatterns.pri \
		../../../../Qt5.7.0/5.7/gcc_64/mkspecs/modules/qt_lib_xmlpatterns_private.pri \
		../../../../Qt5.7.0/5.7/gcc_64/mkspecs/features/qt_functions.prf \
		../../../../Qt5.7.0/5.7/gcc_64/mkspecs/features/qt_config.prf \
		../../../../Qt5.7.0/5.7/gcc_64/mkspecs/linux-g++/qmake.conf \
		../../../../Qt5.7.0/5.7/gcc_64/mkspecs/features/spec_post.prf \
		../Astro/.qmake.stash \
		../../../../Qt5.7.0/5.7/gcc_64/mkspecs/features/exclusive_builds.prf \
		../../../../Qt5.7.0/5.7/gcc_64/mkspecs/features/default_pre.prf \
		../../../../Qt5.7.0/5.7/gcc_64/mkspecs/features/resolve_config.prf \
		../../../../Qt5.7.0/5.7/gcc_64/mkspecs/features/default_post.prf \
		../../../../Qt5.7.0/5.7/gcc_64/mkspecs/features/warn_on.prf \
		../../../../Qt5.7.0/5.7/gcc_64/mkspecs/features/file_copies.prf \
		../../../../Qt5.7.0/5.7/gcc_64/mkspecs/features/testcase_targets.prf \
		../../../../Qt5.7.0/5.7/gcc_64/mkspecs/features/exceptions.prf \
		../../../../Qt5.7.0/5.7/gcc_64/mkspecs/features/yacc.prf \
		../../../../Qt5.7.0/5.7/gcc_64/mkspecs/features/lex.prf \
		../Astro/Astro.pro CommandManager.h \
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
		Distributions.h ../Astro/main.cpp \
		../Astro/CommandManager.cpp \
		../Astro/Models.cpp \
		../Astro/CortexSimulation.cpp \
		../Astro/CortexMeasure.cpp \
		../Astro/Parameters.cpp \
		../Astro/MatrixInverse.cpp \
		../Astro/matrixCholesky.cpp \
		../Astro/LevenbergMarquardt.cpp \
		../Astro/CortexLikelihood.cpp \
		../Astro/Splines.cpp \
		../Astro/Evidence.cpp \
		../Astro/Matrix.cpp \
		../Astro/Distributions.cpp
QMAKE_TARGET  = Astro
DESTDIR       =
TARGET        = Astro


first: all
####### Build rules

$(TARGET):  $(OBJECTS)
	$(LINK) $(LFLAGS) -o $(TARGET) $(OBJECTS) $(OBJCOMP) $(LIBS)

Makefile: ../Astro/Astro.pro ../../../../Qt5.7.0/5.7/gcc_64/mkspecs/linux-g++/qmake.conf ../../../../Qt5.7.0/5.7/gcc_64/mkspecs/features/spec_pre.prf \
../../../../Qt5.7.0/5.7/gcc_64/mkspecs/common/unix.conf \
../../../../Qt5.7.0/5.7/gcc_64/mkspecs/common/linux.conf \
		../../../../Qt5.7.0/5.7/gcc_64/mkspecs/common/sanitize.conf \
		../../../../Qt5.7.0/5.7/gcc_64/mkspecs/common/gcc-base.conf \
		../../../../Qt5.7.0/5.7/gcc_64/mkspecs/common/gcc-base-unix.conf \
		../../../../Qt5.7.0/5.7/gcc_64/mkspecs/common/g++-base.conf \
		../../../../Qt5.7.0/5.7/gcc_64/mkspecs/common/g++-unix.conf \
		../../../../Qt5.7.0/5.7/gcc_64/mkspecs/qconfig.pri \
		../../../../Qt5.7.0/5.7/gcc_64/mkspecs/modules/qt_lib_3dcore.pri \
		../../../../Qt5.7.0/5.7/gcc_64/mkspecs/modules/qt_lib_3dcore_private.pri \
		../../../../Qt5.7.0/5.7/gcc_64/mkspecs/modules/qt_lib_3dextras.pri \
		../../../../Qt5.7.0/5.7/gcc_64/mkspecs/modules/qt_lib_3dextras_private.pri \
		../../../../Qt5.7.0/5.7/gcc_64/mkspecs/modules/qt_lib_3dinput.pri \
		../../../../Qt5.7.0/5.7/gcc_64/mkspecs/modules/qt_lib_3dinput_private.pri \
		../../../../Qt5.7.0/5.7/gcc_64/mkspecs/modules/qt_lib_3dlogic.pri \
		../../../../Qt5.7.0/5.7/gcc_64/mkspecs/modules/qt_lib_3dlogic_private.pri \
		../../../../Qt5.7.0/5.7/gcc_64/mkspecs/modules/qt_lib_3dquick.pri \
		../../../../Qt5.7.0/5.7/gcc_64/mkspecs/modules/qt_lib_3dquick_private.pri \
		../../../../Qt5.7.0/5.7/gcc_64/mkspecs/modules/qt_lib_3dquickextras.pri \
		../../../../Qt5.7.0/5.7/gcc_64/mkspecs/modules/qt_lib_3dquickextras_private.pri \
		../../../../Qt5.7.0/5.7/gcc_64/mkspecs/modules/qt_lib_3dquickinput.pri \
		../../../../Qt5.7.0/5.7/gcc_64/mkspecs/modules/qt_lib_3dquickinput_private.pri \
		../../../../Qt5.7.0/5.7/gcc_64/mkspecs/modules/qt_lib_3dquickrender.pri \
		../../../../Qt5.7.0/5.7/gcc_64/mkspecs/modules/qt_lib_3dquickrender_private.pri \
		../../../../Qt5.7.0/5.7/gcc_64/mkspecs/modules/qt_lib_3drender.pri \
		../../../../Qt5.7.0/5.7/gcc_64/mkspecs/modules/qt_lib_3drender_private.pri \
		../../../../Qt5.7.0/5.7/gcc_64/mkspecs/modules/qt_lib_bluetooth.pri \
		../../../../Qt5.7.0/5.7/gcc_64/mkspecs/modules/qt_lib_bluetooth_private.pri \
		../../../../Qt5.7.0/5.7/gcc_64/mkspecs/modules/qt_lib_bootstrap_private.pri \
		../../../../Qt5.7.0/5.7/gcc_64/mkspecs/modules/qt_lib_charts.pri \
		../../../../Qt5.7.0/5.7/gcc_64/mkspecs/modules/qt_lib_charts_private.pri \
		../../../../Qt5.7.0/5.7/gcc_64/mkspecs/modules/qt_lib_clucene_private.pri \
		../../../../Qt5.7.0/5.7/gcc_64/mkspecs/modules/qt_lib_concurrent.pri \
		../../../../Qt5.7.0/5.7/gcc_64/mkspecs/modules/qt_lib_concurrent_private.pri \
		../../../../Qt5.7.0/5.7/gcc_64/mkspecs/modules/qt_lib_core.pri \
		../../../../Qt5.7.0/5.7/gcc_64/mkspecs/modules/qt_lib_core_private.pri \
		../../../../Qt5.7.0/5.7/gcc_64/mkspecs/modules/qt_lib_datavisualization.pri \
		../../../../Qt5.7.0/5.7/gcc_64/mkspecs/modules/qt_lib_datavisualization_private.pri \
		../../../../Qt5.7.0/5.7/gcc_64/mkspecs/modules/qt_lib_dbus.pri \
		../../../../Qt5.7.0/5.7/gcc_64/mkspecs/modules/qt_lib_dbus_private.pri \
		../../../../Qt5.7.0/5.7/gcc_64/mkspecs/modules/qt_lib_designer.pri \
		../../../../Qt5.7.0/5.7/gcc_64/mkspecs/modules/qt_lib_designer_private.pri \
		../../../../Qt5.7.0/5.7/gcc_64/mkspecs/modules/qt_lib_designercomponents_private.pri \
		../../../../Qt5.7.0/5.7/gcc_64/mkspecs/modules/qt_lib_eglfs_device_lib_private.pri \
		../../../../Qt5.7.0/5.7/gcc_64/mkspecs/modules/qt_lib_gamepad.pri \
		../../../../Qt5.7.0/5.7/gcc_64/mkspecs/modules/qt_lib_gamepad_private.pri \
		../../../../Qt5.7.0/5.7/gcc_64/mkspecs/modules/qt_lib_gui.pri \
		../../../../Qt5.7.0/5.7/gcc_64/mkspecs/modules/qt_lib_gui_private.pri \
		../../../../Qt5.7.0/5.7/gcc_64/mkspecs/modules/qt_lib_help.pri \
		../../../../Qt5.7.0/5.7/gcc_64/mkspecs/modules/qt_lib_help_private.pri \
		../../../../Qt5.7.0/5.7/gcc_64/mkspecs/modules/qt_lib_location.pri \
		../../../../Qt5.7.0/5.7/gcc_64/mkspecs/modules/qt_lib_location_private.pri \
		../../../../Qt5.7.0/5.7/gcc_64/mkspecs/modules/qt_lib_multimedia.pri \
		../../../../Qt5.7.0/5.7/gcc_64/mkspecs/modules/qt_lib_multimedia_private.pri \
		../../../../Qt5.7.0/5.7/gcc_64/mkspecs/modules/qt_lib_multimediawidgets.pri \
		../../../../Qt5.7.0/5.7/gcc_64/mkspecs/modules/qt_lib_multimediawidgets_private.pri \
		../../../../Qt5.7.0/5.7/gcc_64/mkspecs/modules/qt_lib_network.pri \
		../../../../Qt5.7.0/5.7/gcc_64/mkspecs/modules/qt_lib_network_private.pri \
		../../../../Qt5.7.0/5.7/gcc_64/mkspecs/modules/qt_lib_nfc.pri \
		../../../../Qt5.7.0/5.7/gcc_64/mkspecs/modules/qt_lib_nfc_private.pri \
		../../../../Qt5.7.0/5.7/gcc_64/mkspecs/modules/qt_lib_opengl.pri \
		../../../../Qt5.7.0/5.7/gcc_64/mkspecs/modules/qt_lib_opengl_private.pri \
		../../../../Qt5.7.0/5.7/gcc_64/mkspecs/modules/qt_lib_openglextensions.pri \
		../../../../Qt5.7.0/5.7/gcc_64/mkspecs/modules/qt_lib_openglextensions_private.pri \
		../../../../Qt5.7.0/5.7/gcc_64/mkspecs/modules/qt_lib_packetprotocol_private.pri \
		../../../../Qt5.7.0/5.7/gcc_64/mkspecs/modules/qt_lib_platformsupport_private.pri \
		../../../../Qt5.7.0/5.7/gcc_64/mkspecs/modules/qt_lib_positioning.pri \
		../../../../Qt5.7.0/5.7/gcc_64/mkspecs/modules/qt_lib_positioning_private.pri \
		../../../../Qt5.7.0/5.7/gcc_64/mkspecs/modules/qt_lib_printsupport.pri \
		../../../../Qt5.7.0/5.7/gcc_64/mkspecs/modules/qt_lib_printsupport_private.pri \
		../../../../Qt5.7.0/5.7/gcc_64/mkspecs/modules/qt_lib_purchasing.pri \
		../../../../Qt5.7.0/5.7/gcc_64/mkspecs/modules/qt_lib_purchasing_private.pri \
		../../../../Qt5.7.0/5.7/gcc_64/mkspecs/modules/qt_lib_qml.pri \
		../../../../Qt5.7.0/5.7/gcc_64/mkspecs/modules/qt_lib_qml_private.pri \
		../../../../Qt5.7.0/5.7/gcc_64/mkspecs/modules/qt_lib_qmldebug_private.pri \
		../../../../Qt5.7.0/5.7/gcc_64/mkspecs/modules/qt_lib_qmldevtools_private.pri \
		../../../../Qt5.7.0/5.7/gcc_64/mkspecs/modules/qt_lib_qmltest.pri \
		../../../../Qt5.7.0/5.7/gcc_64/mkspecs/modules/qt_lib_qmltest_private.pri \
		../../../../Qt5.7.0/5.7/gcc_64/mkspecs/modules/qt_lib_qtmultimediaquicktools_private.pri \
		../../../../Qt5.7.0/5.7/gcc_64/mkspecs/modules/qt_lib_quick.pri \
		../../../../Qt5.7.0/5.7/gcc_64/mkspecs/modules/qt_lib_quick_private.pri \
		../../../../Qt5.7.0/5.7/gcc_64/mkspecs/modules/qt_lib_quickcontrols2.pri \
		../../../../Qt5.7.0/5.7/gcc_64/mkspecs/modules/qt_lib_quickcontrols2_private.pri \
		../../../../Qt5.7.0/5.7/gcc_64/mkspecs/modules/qt_lib_quickparticles_private.pri \
		../../../../Qt5.7.0/5.7/gcc_64/mkspecs/modules/qt_lib_quicktemplates2_private.pri \
		../../../../Qt5.7.0/5.7/gcc_64/mkspecs/modules/qt_lib_quickwidgets.pri \
		../../../../Qt5.7.0/5.7/gcc_64/mkspecs/modules/qt_lib_quickwidgets_private.pri \
		../../../../Qt5.7.0/5.7/gcc_64/mkspecs/modules/qt_lib_script.pri \
		../../../../Qt5.7.0/5.7/gcc_64/mkspecs/modules/qt_lib_script_private.pri \
		../../../../Qt5.7.0/5.7/gcc_64/mkspecs/modules/qt_lib_scripttools.pri \
		../../../../Qt5.7.0/5.7/gcc_64/mkspecs/modules/qt_lib_scripttools_private.pri \
		../../../../Qt5.7.0/5.7/gcc_64/mkspecs/modules/qt_lib_scxml.pri \
		../../../../Qt5.7.0/5.7/gcc_64/mkspecs/modules/qt_lib_scxml_private.pri \
		../../../../Qt5.7.0/5.7/gcc_64/mkspecs/modules/qt_lib_sensors.pri \
		../../../../Qt5.7.0/5.7/gcc_64/mkspecs/modules/qt_lib_sensors_private.pri \
		../../../../Qt5.7.0/5.7/gcc_64/mkspecs/modules/qt_lib_serialbus.pri \
		../../../../Qt5.7.0/5.7/gcc_64/mkspecs/modules/qt_lib_serialbus_private.pri \
		../../../../Qt5.7.0/5.7/gcc_64/mkspecs/modules/qt_lib_serialport.pri \
		../../../../Qt5.7.0/5.7/gcc_64/mkspecs/modules/qt_lib_serialport_private.pri \
		../../../../Qt5.7.0/5.7/gcc_64/mkspecs/modules/qt_lib_sql.pri \
		../../../../Qt5.7.0/5.7/gcc_64/mkspecs/modules/qt_lib_sql_private.pri \
		../../../../Qt5.7.0/5.7/gcc_64/mkspecs/modules/qt_lib_svg.pri \
		../../../../Qt5.7.0/5.7/gcc_64/mkspecs/modules/qt_lib_svg_private.pri \
		../../../../Qt5.7.0/5.7/gcc_64/mkspecs/modules/qt_lib_testlib.pri \
		../../../../Qt5.7.0/5.7/gcc_64/mkspecs/modules/qt_lib_testlib_private.pri \
		../../../../Qt5.7.0/5.7/gcc_64/mkspecs/modules/qt_lib_uiplugin.pri \
		../../../../Qt5.7.0/5.7/gcc_64/mkspecs/modules/qt_lib_uitools.pri \
		../../../../Qt5.7.0/5.7/gcc_64/mkspecs/modules/qt_lib_uitools_private.pri \
		../../../../Qt5.7.0/5.7/gcc_64/mkspecs/modules/qt_lib_webchannel.pri \
		../../../../Qt5.7.0/5.7/gcc_64/mkspecs/modules/qt_lib_webchannel_private.pri \
		../../../../Qt5.7.0/5.7/gcc_64/mkspecs/modules/qt_lib_webengine.pri \
		../../../../Qt5.7.0/5.7/gcc_64/mkspecs/modules/qt_lib_webengine_private.pri \
		../../../../Qt5.7.0/5.7/gcc_64/mkspecs/modules/qt_lib_webenginecore.pri \
		../../../../Qt5.7.0/5.7/gcc_64/mkspecs/modules/qt_lib_webenginecore_private.pri \
		../../../../Qt5.7.0/5.7/gcc_64/mkspecs/modules/qt_lib_webenginecoreheaders_private.pri \
		../../../../Qt5.7.0/5.7/gcc_64/mkspecs/modules/qt_lib_webenginewidgets.pri \
		../../../../Qt5.7.0/5.7/gcc_64/mkspecs/modules/qt_lib_webenginewidgets_private.pri \
		../../../../Qt5.7.0/5.7/gcc_64/mkspecs/modules/qt_lib_websockets.pri \
		../../../../Qt5.7.0/5.7/gcc_64/mkspecs/modules/qt_lib_websockets_private.pri \
		../../../../Qt5.7.0/5.7/gcc_64/mkspecs/modules/qt_lib_webview.pri \
		../../../../Qt5.7.0/5.7/gcc_64/mkspecs/modules/qt_lib_webview_private.pri \
		../../../../Qt5.7.0/5.7/gcc_64/mkspecs/modules/qt_lib_widgets.pri \
		../../../../Qt5.7.0/5.7/gcc_64/mkspecs/modules/qt_lib_widgets_private.pri \
		../../../../Qt5.7.0/5.7/gcc_64/mkspecs/modules/qt_lib_x11extras.pri \
		../../../../Qt5.7.0/5.7/gcc_64/mkspecs/modules/qt_lib_x11extras_private.pri \
		../../../../Qt5.7.0/5.7/gcc_64/mkspecs/modules/qt_lib_xcb_qpa_lib_private.pri \
		../../../../Qt5.7.0/5.7/gcc_64/mkspecs/modules/qt_lib_xml.pri \
		../../../../Qt5.7.0/5.7/gcc_64/mkspecs/modules/qt_lib_xml_private.pri \
		../../../../Qt5.7.0/5.7/gcc_64/mkspecs/modules/qt_lib_xmlpatterns.pri \
		../../../../Qt5.7.0/5.7/gcc_64/mkspecs/modules/qt_lib_xmlpatterns_private.pri \
		../../../../Qt5.7.0/5.7/gcc_64/mkspecs/features/qt_functions.prf \
		../../../../Qt5.7.0/5.7/gcc_64/mkspecs/features/qt_config.prf \
		../../../../Qt5.7.0/5.7/gcc_64/mkspecs/linux-g++/qmake.conf \
		../../../../Qt5.7.0/5.7/gcc_64/mkspecs/features/spec_post.prf \
		.qmake.stash \
		../../../../Qt5.7.0/5.7/gcc_64/mkspecs/features/exclusive_builds.prf \
		../../../../Qt5.7.0/5.7/gcc_64/mkspecs/features/default_pre.prf \
		../../../../Qt5.7.0/5.7/gcc_64/mkspecs/features/resolve_config.prf \
		../../../../Qt5.7.0/5.7/gcc_64/mkspecs/features/default_post.prf \
		../../../../Qt5.7.0/5.7/gcc_64/mkspecs/features/warn_on.prf \
		../../../../Qt5.7.0/5.7/gcc_64/mkspecs/features/file_copies.prf \
		../../../../Qt5.7.0/5.7/gcc_64/mkspecs/features/testcase_targets.prf \
		../../../../Qt5.7.0/5.7/gcc_64/mkspecs/features/exceptions.prf \
		../../../../Qt5.7.0/5.7/gcc_64/mkspecs/features/yacc.prf \
		../../../../Qt5.7.0/5.7/gcc_64/mkspecs/features/lex.prf \
		../Astro/Astro.pro
	$(QMAKE) -spec linux-g++ -o Makefile ../Astro/Astro.pro


all: Makefile $(TARGET)

dist: distdir FORCE
	(cd `dirname $(DISTDIR)` && $(TAR) $(DISTNAME).tar $(DISTNAME) && $(COMPRESS) $(DISTNAME).tar) && $(MOVE) `dirname $(DISTDIR)`/$(DISTNAME).tar.gz . && $(DEL_FILE) -r $(DISTDIR)

distdir: FORCE
@test -d $(DISTDIR) || mkdir -p $(DISTDIR)
$(COPY_FILE) --parents $(DIST) $(DISTDIR)/


clean: compiler_clean
-$(DEL_FILE) $(OBJECTS)
-$(DEL_FILE) *~ core *.core


distclean: clean
-$(DEL_FILE) $(TARGET)
-$(DEL_FILE) .qmake.stash
	-$(DEL_FILE) Makefile


####### Sub-libraries

check: first

benchmark: first

compiler_yacc_decl_make_all:
compiler_yacc_decl_clean:
compiler_yacc_impl_make_all:
compiler_yacc_impl_clean:
compiler_lex_make_all:
compiler_lex_clean:
compiler_clean:

####### Compile

main.o: ../Astro/main.cpp ../Astro/CommandManager.h \
../Astro/BaseClass.h \
../Astro/CortexMeasure.h \
		../Astro/LevenbergMarquardt.h \
		../Astro/Parameters.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o main.o ../Astro/main.cpp

CommandManager.o: ../Astro/CommandManager.cpp ../Astro/Models.h \
../Astro/CortexMeasure.h \
../Astro/LevenbergMarquardt.h \
		../Astro/Parameters.h \
		../Astro/BaseClass.h \
		../Astro/CortexSimulation.h \
		../Astro/CommandManager.h \
		../Astro/CortexLikelihood.h \
		../Astro/Evidence.h \
		../Astro/Matrix.h \
		../Astro/Distributions.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o CommandManager.o ../Astro/CommandManager.cpp

Models.o: ../Astro/Models.cpp ../Astro/Models.h \
../Astro/CortexMeasure.h \
../Astro/LevenbergMarquardt.h \
		../Astro/Parameters.h \
		../Astro/BaseClass.h \
		../Astro/CortexSimulation.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o Models.o ../Astro/Models.cpp

CortexSimulation.o: ../Astro/CortexSimulation.cpp ../Astro/CortexSimulation.h \
../Astro/Parameters.h \
../Astro/BaseClass.h \
		../Astro/CommandManager.h \
		../Astro/CortexMeasure.h \
		../Astro/LevenbergMarquardt.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o CortexSimulation.o ../Astro/CortexSimulation.cpp

CortexMeasure.o: ../Astro/CortexMeasure.cpp ../Astro/CommandManager.h \
../Astro/BaseClass.h \
../Astro/CortexMeasure.h \
		../Astro/LevenbergMarquardt.h \
		../Astro/Parameters.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o CortexMeasure.o ../Astro/CortexMeasure.cpp

Parameters.o: ../Astro/Parameters.cpp ../Astro/Parameters.h \
../Astro/BaseClass.h \
../Astro/MatrixInverse.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o Parameters.o ../Astro/Parameters.cpp

MatrixInverse.o: ../Astro/MatrixInverse.cpp ../Astro/MatrixInverse.h
$(CXX) -c $(CXXFLAGS) $(INCPATH) -o MatrixInverse.o ../Astro/MatrixInverse.cpp

matrixCholesky.o: ../Astro/matrixCholesky.cpp ../Astro/MatrixInverse.h
$(CXX) -c $(CXXFLAGS) $(INCPATH) -o matrixCholesky.o ../Astro/matrixCholesky.cpp

LevenbergMarquardt.o: ../Astro/LevenbergMarquardt.cpp ../Astro/LevenbergMarquardt.h \
../Astro/Parameters.h \
../Astro/BaseClass.h \
		../Astro/MatrixInverse.h \
		../Astro/CortexLikelihood.h \
		../Astro/Models.h \
		../Astro/CortexMeasure.h \
		../Astro/CortexSimulation.h \
		../Astro/Evidence.h \
		../Astro/Matrix.h \
		../Astro/Distributions.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o LevenbergMarquardt.o ../Astro/LevenbergMarquardt.cpp

CortexLikelihood.o: ../Astro/CortexLikelihood.cpp ../Astro/CortexLikelihood.h \
../Astro/Models.h \
../Astro/CortexMeasure.h \
		../Astro/LevenbergMarquardt.h \
		../Astro/Parameters.h \
		../Astro/BaseClass.h \
		../Astro/CortexSimulation.h \
		../Astro/Evidence.h \
		../Astro/Matrix.h \
		../Astro/Distributions.h \
		../Astro/CommandManager.h \
		../Astro/Splines.h \
		../Astro/MatrixInverse.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o CortexLikelihood.o ../Astro/CortexLikelihood.cpp

Splines.o: ../Astro/Splines.cpp ../Astro/Splines.h \
../Astro/MatrixInverse.h
$(CXX) -c $(CXXFLAGS) $(INCPATH) -o Splines.o ../Astro/Splines.cpp

Evidence.o: ../Astro/Evidence.cpp ../Astro/Evidence.h \
../Astro/Matrix.h \
../Astro/Distributions.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o Evidence.o ../Astro/Evidence.cpp

Matrix.o: ../Astro/Matrix.cpp ../Astro/Matrix.h
$(CXX) -c $(CXXFLAGS) $(INCPATH) -o Matrix.o ../Astro/Matrix.cpp

Distributions.o: ../Astro/Distributions.cpp ../Astro/Distributions.h \
../Astro/Matrix.h
$(CXX) -c $(CXXFLAGS) $(INCPATH) -o Distributions.o ../Astro/Distributions.cpp

####### Install

install:  FORCE

uninstall:  FORCE

FORCE:

