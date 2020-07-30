#!/bin/bash

g++ exKF2.cpp -o exKF2 -llapack -lblas -lm 
./exKF2 > exKF2.dat

gmt psxy -R0/1500/0/2000 -JX20/10 -Bnesw -Bx100 -By200 -T -K > two.ps

cat exKF2.dat|awk '{print $3,$4}'|\
gmt psxy -R -J -B -Sc0.02c -Gyellow -K -O >> two.ps;
cat exKF2.dat|awk '{print $5,$6}'|\
gmt psxy -R -J -B -Sc0.02c -Gred -K -O >> two.ps;
cat exKF2.dat|awk '{print $1,$2}'|\
gmt psxy -R -J -B -Sc0.02c -Gblue -K -O >> two.ps;
gmt psxy -R -J -B -T -O >> two.ps;
gmt psconvert two.ps -A -Tg -E600;