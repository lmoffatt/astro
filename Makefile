#############################################################################
# Makefile for building: Astro
# Generated by qmake (3.0) (Qt 5.5.0)
# Project:  ../Astro/Astro.pro
# Template: app
# Command: /home/luciano/Qt5.5.0/5.5/gcc_64/bin/qmake -spec linux-g++ -o Makefile ../Astro/Astro.pro
#############################################################################

MAKEFILE      = Makefile

####### Compiler, tools and options

CC            = gcc
CXX           = g++-5
DEFINES       = 
CFLAGS        = -pipe -O2 -Wall -W -fPIC $(DEFINES)
CXXFLAGS      = -pipe -O2 -std=c++11 -Wall -W -fPIC $(DEFINES)
INCPATH       = -I../Astro -I. -I../../../../Qt5.5.0/5.5/gcc_64/mkspecs/linux-g++
QMAKE         = /home/luciano/Qt5.5.0/5.5/gcc_64/bin/qmake
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
DISTDIR = /home/luciano/Data/celulas/Astro/build-Astro-Desktop-Release/.tmp/Astro1.0.0
LINK          = g++
LFLAGS        = -Wl,-O1 -Wl,-rpath,/home/luciano/Qt5.5.0/5.5/gcc_64
LIBS          = $(SUBLIBS) -L/home/luciano/Data/celulas/Astro/Astro/bin -lblas -llapack 
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
		../Astro/BayesIteration.cpp \
		../Astro/CortexLikelihood.cpp 
OBJECTS       = main.o \
		CommandManager.o \
		Models.o \
		CortexSimulation.o \
		CortexMeasure.o \
		Parameters.o \
		MatrixInverse.o \
		matrixCholesky.o \
		LevenbergMarquardt.o \
		BayesIteration.o \
		CortexLikelihood.o
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
		../Astro/run/parameters_011.txt \
		../Astro/run/parameters_012.txt \
		../Astro/run/parameters_013.txt \
		../Astro/run/.txt \
		../Astro/run/parameters_013.txt \
		../Astro/run/parameters_111.txt \
		../Astro/run/parameters_112.txt \
		../Astro/run/parameters_113.txt \
		../Astro/run/parameters_114.txt \
		../Astro/run/parameters_116.txt \
		../Astro/Makefile \
		../../../../Qt5.5.0/5.5/gcc_64/mkspecs/features/spec_pre.prf \
		../../../../Qt5.5.0/5.5/gcc_64/mkspecs/common/unix.conf \
		../../../../Qt5.5.0/5.5/gcc_64/mkspecs/common/linux.conf \
		../../../../Qt5.5.0/5.5/gcc_64/mkspecs/common/sanitize.conf \
		../../../../Qt5.5.0/5.5/gcc_64/mkspecs/common/gcc-base.conf \
		../../../../Qt5.5.0/5.5/gcc_64/mkspecs/common/gcc-base-unix.conf \
		../../../../Qt5.5.0/5.5/gcc_64/mkspecs/common/g++-base.conf \
		../../../../Qt5.5.0/5.5/gcc_64/mkspecs/common/g++-unix.conf \
		../../../../Qt5.5.0/5.5/gcc_64/mkspecs/qconfig.pri \
		../../../../Qt5.5.0/5.5/gcc_64/mkspecs/modules/qt_lib_3dcore.pri \
		../../../../Qt5.5.0/5.5/gcc_64/mkspecs/modules/qt_lib_3dcore_private.pri \
		../../../../Qt5.5.0/5.5/gcc_64/mkspecs/modules/qt_lib_3dinput.pri \
		../../../../Qt5.5.0/5.5/gcc_64/mkspecs/modules/qt_lib_3dinput_private.pri \
		../../../../Qt5.5.0/5.5/gcc_64/mkspecs/modules/qt_lib_3dquick.pri \
		../../../../Qt5.5.0/5.5/gcc_64/mkspecs/modules/qt_lib_3dquick_private.pri \
		../../../../Qt5.5.0/5.5/gcc_64/mkspecs/modules/qt_lib_3dquickrenderer.pri \
		../../../../Qt5.5.0/5.5/gcc_64/mkspecs/modules/qt_lib_3dquickrenderer_private.pri \
		../../../../Qt5.5.0/5.5/gcc_64/mkspecs/modules/qt_lib_3drenderer.pri \
		../../../../Qt5.5.0/5.5/gcc_64/mkspecs/modules/qt_lib_3drenderer_private.pri \
		../../../../Qt5.5.0/5.5/gcc_64/mkspecs/modules/qt_lib_bluetooth.pri \
		../../../../Qt5.5.0/5.5/gcc_64/mkspecs/modules/qt_lib_bluetooth_private.pri \
		../../../../Qt5.5.0/5.5/gcc_64/mkspecs/modules/qt_lib_bootstrap_private.pri \
		../../../../Qt5.5.0/5.5/gcc_64/mkspecs/modules/qt_lib_clucene_private.pri \
		../../../../Qt5.5.0/5.5/gcc_64/mkspecs/modules/qt_lib_concurrent.pri \
		../../../../Qt5.5.0/5.5/gcc_64/mkspecs/modules/qt_lib_concurrent_private.pri \
		../../../../Qt5.5.0/5.5/gcc_64/mkspecs/modules/qt_lib_core.pri \
		../../../../Qt5.5.0/5.5/gcc_64/mkspecs/modules/qt_lib_core_private.pri \
		../../../../Qt5.5.0/5.5/gcc_64/mkspecs/modules/qt_lib_dbus.pri \
		../../../../Qt5.5.0/5.5/gcc_64/mkspecs/modules/qt_lib_dbus_private.pri \
		../../../../Qt5.5.0/5.5/gcc_64/mkspecs/modules/qt_lib_declarative.pri \
		../../../../Qt5.5.0/5.5/gcc_64/mkspecs/modules/qt_lib_declarative_private.pri \
		../../../../Qt5.5.0/5.5/gcc_64/mkspecs/modules/qt_lib_designer.pri \
		../../../../Qt5.5.0/5.5/gcc_64/mkspecs/modules/qt_lib_designer_private.pri \
		../../../../Qt5.5.0/5.5/gcc_64/mkspecs/modules/qt_lib_designercomponents_private.pri \
		../../../../Qt5.5.0/5.5/gcc_64/mkspecs/modules/qt_lib_enginio.pri \
		../../../../Qt5.5.0/5.5/gcc_64/mkspecs/modules/qt_lib_enginio_private.pri \
		../../../../Qt5.5.0/5.5/gcc_64/mkspecs/modules/qt_lib_gui.pri \
		../../../../Qt5.5.0/5.5/gcc_64/mkspecs/modules/qt_lib_gui_private.pri \
		../../../../Qt5.5.0/5.5/gcc_64/mkspecs/modules/qt_lib_help.pri \
		../../../../Qt5.5.0/5.5/gcc_64/mkspecs/modules/qt_lib_help_private.pri \
		../../../../Qt5.5.0/5.5/gcc_64/mkspecs/modules/qt_lib_location.pri \
		../../../../Qt5.5.0/5.5/gcc_64/mkspecs/modules/qt_lib_location_private.pri \
		../../../../Qt5.5.0/5.5/gcc_64/mkspecs/modules/qt_lib_multimedia.pri \
		../../../../Qt5.5.0/5.5/gcc_64/mkspecs/modules/qt_lib_multimedia_private.pri \
		../../../../Qt5.5.0/5.5/gcc_64/mkspecs/modules/qt_lib_multimediawidgets.pri \
		../../../../Qt5.5.0/5.5/gcc_64/mkspecs/modules/qt_lib_multimediawidgets_private.pri \
		../../../../Qt5.5.0/5.5/gcc_64/mkspecs/modules/qt_lib_network.pri \
		../../../../Qt5.5.0/5.5/gcc_64/mkspecs/modules/qt_lib_network_private.pri \
		../../../../Qt5.5.0/5.5/gcc_64/mkspecs/modules/qt_lib_nfc.pri \
		../../../../Qt5.5.0/5.5/gcc_64/mkspecs/modules/qt_lib_nfc_private.pri \
		../../../../Qt5.5.0/5.5/gcc_64/mkspecs/modules/qt_lib_opengl.pri \
		../../../../Qt5.5.0/5.5/gcc_64/mkspecs/modules/qt_lib_opengl_private.pri \
		../../../../Qt5.5.0/5.5/gcc_64/mkspecs/modules/qt_lib_openglextensions.pri \
		../../../../Qt5.5.0/5.5/gcc_64/mkspecs/modules/qt_lib_openglextensions_private.pri \
		../../../../Qt5.5.0/5.5/gcc_64/mkspecs/modules/qt_lib_platformsupport_private.pri \
		../../../../Qt5.5.0/5.5/gcc_64/mkspecs/modules/qt_lib_positioning.pri \
		../../../../Qt5.5.0/5.5/gcc_64/mkspecs/modules/qt_lib_positioning_private.pri \
		../../../../Qt5.5.0/5.5/gcc_64/mkspecs/modules/qt_lib_printsupport.pri \
		../../../../Qt5.5.0/5.5/gcc_64/mkspecs/modules/qt_lib_printsupport_private.pri \
		../../../../Qt5.5.0/5.5/gcc_64/mkspecs/modules/qt_lib_qml.pri \
		../../../../Qt5.5.0/5.5/gcc_64/mkspecs/modules/qt_lib_qml_private.pri \
		../../../../Qt5.5.0/5.5/gcc_64/mkspecs/modules/qt_lib_qmldevtools_private.pri \
		../../../../Qt5.5.0/5.5/gcc_64/mkspecs/modules/qt_lib_qmltest.pri \
		../../../../Qt5.5.0/5.5/gcc_64/mkspecs/modules/qt_lib_qmltest_private.pri \
		../../../../Qt5.5.0/5.5/gcc_64/mkspecs/modules/qt_lib_qtmultimediaquicktools_private.pri \
		../../../../Qt5.5.0/5.5/gcc_64/mkspecs/modules/qt_lib_quick.pri \
		../../../../Qt5.5.0/5.5/gcc_64/mkspecs/modules/qt_lib_quick_private.pri \
		../../../../Qt5.5.0/5.5/gcc_64/mkspecs/modules/qt_lib_quickparticles_private.pri \
		../../../../Qt5.5.0/5.5/gcc_64/mkspecs/modules/qt_lib_quickwidgets.pri \
		../../../../Qt5.5.0/5.5/gcc_64/mkspecs/modules/qt_lib_quickwidgets_private.pri \
		../../../../Qt5.5.0/5.5/gcc_64/mkspecs/modules/qt_lib_script.pri \
		../../../../Qt5.5.0/5.5/gcc_64/mkspecs/modules/qt_lib_script_private.pri \
		../../../../Qt5.5.0/5.5/gcc_64/mkspecs/modules/qt_lib_scripttools.pri \
		../../../../Qt5.5.0/5.5/gcc_64/mkspecs/modules/qt_lib_scripttools_private.pri \
		../../../../Qt5.5.0/5.5/gcc_64/mkspecs/modules/qt_lib_sensors.pri \
		../../../../Qt5.5.0/5.5/gcc_64/mkspecs/modules/qt_lib_sensors_private.pri \
		../../../../Qt5.5.0/5.5/gcc_64/mkspecs/modules/qt_lib_serialport.pri \
		../../../../Qt5.5.0/5.5/gcc_64/mkspecs/modules/qt_lib_serialport_private.pri \
		../../../../Qt5.5.0/5.5/gcc_64/mkspecs/modules/qt_lib_sql.pri \
		../../../../Qt5.5.0/5.5/gcc_64/mkspecs/modules/qt_lib_sql_private.pri \
		../../../../Qt5.5.0/5.5/gcc_64/mkspecs/modules/qt_lib_svg.pri \
		../../../../Qt5.5.0/5.5/gcc_64/mkspecs/modules/qt_lib_svg_private.pri \
		../../../../Qt5.5.0/5.5/gcc_64/mkspecs/modules/qt_lib_testlib.pri \
		../../../../Qt5.5.0/5.5/gcc_64/mkspecs/modules/qt_lib_testlib_private.pri \
		../../../../Qt5.5.0/5.5/gcc_64/mkspecs/modules/qt_lib_uiplugin.pri \
		../../../../Qt5.5.0/5.5/gcc_64/mkspecs/modules/qt_lib_uitools.pri \
		../../../../Qt5.5.0/5.5/gcc_64/mkspecs/modules/qt_lib_uitools_private.pri \
		../../../../Qt5.5.0/5.5/gcc_64/mkspecs/modules/qt_lib_webchannel.pri \
		../../../../Qt5.5.0/5.5/gcc_64/mkspecs/modules/qt_lib_webchannel_private.pri \
		../../../../Qt5.5.0/5.5/gcc_64/mkspecs/modules/qt_lib_webengine.pri \
		../../../../Qt5.5.0/5.5/gcc_64/mkspecs/modules/qt_lib_webengine_private.pri \
		../../../../Qt5.5.0/5.5/gcc_64/mkspecs/modules/qt_lib_webenginecore.pri \
		../../../../Qt5.5.0/5.5/gcc_64/mkspecs/modules/qt_lib_webenginecore_private.pri \
		../../../../Qt5.5.0/5.5/gcc_64/mkspecs/modules/qt_lib_webenginewidgets.pri \
		../../../../Qt5.5.0/5.5/gcc_64/mkspecs/modules/qt_lib_webenginewidgets_private.pri \
		../../../../Qt5.5.0/5.5/gcc_64/mkspecs/modules/qt_lib_webkit.pri \
		../../../../Qt5.5.0/5.5/gcc_64/mkspecs/modules/qt_lib_webkit_private.pri \
		../../../../Qt5.5.0/5.5/gcc_64/mkspecs/modules/qt_lib_webkitwidgets.pri \
		../../../../Qt5.5.0/5.5/gcc_64/mkspecs/modules/qt_lib_webkitwidgets_private.pri \
		../../../../Qt5.5.0/5.5/gcc_64/mkspecs/modules/qt_lib_websockets.pri \
		../../../../Qt5.5.0/5.5/gcc_64/mkspecs/modules/qt_lib_websockets_private.pri \
		../../../../Qt5.5.0/5.5/gcc_64/mkspecs/modules/qt_lib_webview_private.pri \
		../../../../Qt5.5.0/5.5/gcc_64/mkspecs/modules/qt_lib_widgets.pri \
		../../../../Qt5.5.0/5.5/gcc_64/mkspecs/modules/qt_lib_widgets_private.pri \
		../../../../Qt5.5.0/5.5/gcc_64/mkspecs/modules/qt_lib_x11extras.pri \
		../../../../Qt5.5.0/5.5/gcc_64/mkspecs/modules/qt_lib_x11extras_private.pri \
		../../../../Qt5.5.0/5.5/gcc_64/mkspecs/modules/qt_lib_xcb_qpa_lib_private.pri \
		../../../../Qt5.5.0/5.5/gcc_64/mkspecs/modules/qt_lib_xml.pri \
		../../../../Qt5.5.0/5.5/gcc_64/mkspecs/modules/qt_lib_xml_private.pri \
		../../../../Qt5.5.0/5.5/gcc_64/mkspecs/modules/qt_lib_xmlpatterns.pri \
		../../../../Qt5.5.0/5.5/gcc_64/mkspecs/modules/qt_lib_xmlpatterns_private.pri \
		../../../../Qt5.5.0/5.5/gcc_64/mkspecs/features/qt_functions.prf \
		../../../../Qt5.5.0/5.5/gcc_64/mkspecs/features/qt_config.prf \
		../../../../Qt5.5.0/5.5/gcc_64/mkspecs/linux-g++/qmake.conf \
		../../../../Qt5.5.0/5.5/gcc_64/mkspecs/features/spec_post.prf \
		../../../../Qt5.5.0/5.5/gcc_64/mkspecs/features/exclusive_builds.prf \
		../../../../Qt5.5.0/5.5/gcc_64/mkspecs/features/default_pre.prf \
		../../../../Qt5.5.0/5.5/gcc_64/mkspecs/features/resolve_config.prf \
		../../../../Qt5.5.0/5.5/gcc_64/mkspecs/features/default_post.prf \
		../../../../Qt5.5.0/5.5/gcc_64/mkspecs/features/warn_on.prf \
		../../../../Qt5.5.0/5.5/gcc_64/mkspecs/features/testcase_targets.prf \
		../../../../Qt5.5.0/5.5/gcc_64/mkspecs/features/exceptions.prf \
		../../../../Qt5.5.0/5.5/gcc_64/mkspecs/features/yacc.prf \
		../../../../Qt5.5.0/5.5/gcc_64/mkspecs/features/lex.prf \
		../Astro/Astro.pro CommandManager.h \
		Models.h \
		CortexSimulation.h \
		CortexMeasure.h \
		Parameters.h \
		MatrixInverse.h \
		LevenbergMarquardt.h \
		BayesIteration.h \
		CortexLikelihood.h \
		BaseClass.h \
		MCMC.h ../Astro/main.cpp \
		../Astro/CommandManager.cpp \
		../Astro/Models.cpp \
		../Astro/CortexSimulation.cpp \
		../Astro/CortexMeasure.cpp \
		../Astro/Parameters.cpp \
		../Astro/MatrixInverse.cpp \
		../Astro/matrixCholesky.cpp \
		../Astro/LevenbergMarquardt.cpp \
		../Astro/BayesIteration.cpp \
		../Astro/CortexLikelihood.cpp
QMAKE_TARGET  = Astro
DESTDIR       = #avoid trailing-slash linebreak
TARGET        = Astro


first: all
####### Implicit rules

.SUFFIXES: .o .c .cpp .cc .cxx .C

.cpp.o:
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o "$@" "$<"

.cc.o:
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o "$@" "$<"

.cxx.o:
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o "$@" "$<"

.C.o:
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o "$@" "$<"

.c.o:
	$(CC) -c $(CFLAGS) $(INCPATH) -o "$@" "$<"

####### Build rules

$(TARGET):  $(OBJECTS)  
	$(LINK) $(LFLAGS) -o $(TARGET) $(OBJECTS) $(OBJCOMP) $(LIBS)


qmake_all: FORCE


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
	-$(DEL_FILE) Makefile


####### Sub-libraries

check: first

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
		../Astro/BayesIteration.h \
		../Astro/CommandManager.h \
		../Astro/CortexLikelihood.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o CommandManager.o ../Astro/CommandManager.cpp

Models.o: ../Astro/Models.cpp ../Astro/Models.h \
		../Astro/CortexMeasure.h \
		../Astro/LevenbergMarquardt.h \
		../Astro/Parameters.h \
		../Astro/BaseClass.h \
		../Astro/CortexSimulation.h \
		../Astro/BayesIteration.h
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
		../Astro/BayesIteration.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o LevenbergMarquardt.o ../Astro/LevenbergMarquardt.cpp

BayesIteration.o: ../Astro/BayesIteration.cpp ../Astro/BayesIteration.h \
		../Astro/LevenbergMarquardt.h \
		../Astro/Parameters.h \
		../Astro/BaseClass.h \
		../Astro/MatrixInverse.h \
		../Astro/CortexLikelihood.h \
		../Astro/Models.h \
		../Astro/CortexMeasure.h \
		../Astro/CortexSimulation.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o BayesIteration.o ../Astro/BayesIteration.cpp

CortexLikelihood.o: ../Astro/CortexLikelihood.cpp ../Astro/CortexLikelihood.h \
		../Astro/Models.h \
		../Astro/CortexMeasure.h \
		../Astro/LevenbergMarquardt.h \
		../Astro/Parameters.h \
		../Astro/BaseClass.h \
		../Astro/CortexSimulation.h \
		../Astro/BayesIteration.h \
		../Astro/CommandManager.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o CortexLikelihood.o ../Astro/CortexLikelihood.cpp

####### Install

install:  FORCE

uninstall:  FORCE

FORCE:

