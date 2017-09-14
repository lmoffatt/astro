#!/bin/bash

cd /home/lmoffatt.inquimae/Astro/Astro/run/m01
qsub p_m01
qsub p_m01_1
qsub p_m01_2

cd /home/lmoffatt.inquimae/Astro/Astro/run/m10
qsub p_m10
qsub p_m10_1
qsub p_m10_2

cd /home/lmoffatt.inquimae/Astro/Astro/run/m101
qsub p_m101
qsub p_m101_1
qsub p_m101_2

