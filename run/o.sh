#!/bin/bash
#
# Script para correr trabajo serial
#

# Opciones SGE
#$ -l h_rt=5:00:00 # Setea 1 horas de wall clock time
#$ -V  # Exporta las variables de entorno
#$ -N test  # El nombre del job
#$ -cwd # Cambia al directorio actual

# Comando para correr el programa, tal cual lo llamaríamos desde la línea de comandos
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/home/lmoffatt.inquimae/opt/gcc-5.2.0/lib64/

/home/lmoffatt.inquimae/Astro/release/Astro opt_script LABEL parameters_zero
