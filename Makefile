#############################################################################
# Makefile for building: Astro
# Generated by qmake (3.0) (Qt 5.5.0)
# Project:  ../Astro/Astro.pro
# Template: app
# Command: /home/luciano/Qt5.5.0/5.5/gcc_64/bin/qmake -spec linux-g++ -o Makefile ../Astro/Astro.pro
#############################################################################

MAKEFILE      = Makefile

####### Compiler, tools and options
CC            = ~/opt/gcc-5.2.0/bin/gcc
CXX           = ~/opt/gcc-5.2.0/bin/g++
DEFINES       =
CFLAGS        = -m64 -pipe -O2 -Wall -W $(DEFINES)
CXXFLAGS      = -m64 -pipe -std=c++11 -O2 -Wall -W $(DEFINES)
INCPATH       = -I../Astro -I. -I../../../../Qt5.5.0/5.5/gcc_64/mkspecs/linux-g++
QMAKE         = /home/luciano/Qt5.5.0/5.5/gcc_64/bin/qmake
DEL_FILE      = rm -f

LFLAGS        = -Wl,-O1 -Wl,-rpath, /home/lmoffatt.inquimae/opt/gcc-5.2.0/lib/  /home/lmoffatt.inquimae/opt/gcc-5.2.0/lib64
LINK          = ~/opt/gcc-5.2.0/bin/g++

LIBS          = $(SUBLIBS)   -lblas -llapack

AR            = ar cqs
RANLIB        =
QMAKE         = /home/luciano/qt-everywhere-opensource-src-4.8.3/bin/qmake
TAR           = tar -cf
COMPRESS      = gzip -9f
COPY          = cp -f
SED           = sed
COPY_FILE     = $(COPY)
COPY_DIR      = $(COPY) -r
STRIP         = strip
INSTALL_FILE  = install -m 644 -p
INSTALL_DIR   = $(COPY_DIR)
INSTALL_PROGRAM = install -m 755 -p
DEL_FILE      = rm -f
SYMLINK       = ln -f -s
DEL_DIR       = rmdir
MOVE          = mv -f
CHK_DIR_EXISTS= test -d
MKDIR         = mkdir -p

DISTNAME      = Astro1.0.0
DISTDIR = ~/Astro/release/.tmp/Astro1.0.0
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

