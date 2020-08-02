#!/bin/bash

for i in 0 1 2;do
    g++ exKF${i}.cpp -o exKF${i} -llapack -lblas -lm 
    ./exKF${i} > exKF${i}.dat

    gmt psxy -R0/1500/0/2000 -JX20/10 -Bnesw -Bx100 -By200 -T -K > pic${i}.ps
    cat exKF${i}.dat|awk '{print $3,$4}'|\
        gmt psxy -R -J -B -Sc0.02c -Gyellow -K -O >> pic${i}.ps;
    cat exKF${i}.dat|awk '{print $5,$6}'|\
        gmt psxy -R -J -B -Sc0.02c -Gred -K -O >> pic${i}.ps;
    cat exKF${i}.dat|awk '{print $1,$2}'|\
    gmt psxy -R -J -B -Sc0.02c -Gblue -K -O >> pic${i}.ps;
    gmt psxy -R -J -B -T -O >> pic${i}.ps;
    gmt psconvert pic${i}.ps -A -Tg -E600;

    rm exKF${i} exKF${i}.dat pic${i}.ps

done