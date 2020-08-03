#!/bin/bash

  for i in 0  1 2;do
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
for i in 0 1;do
g++ TDOA${i}.cpp -o TDOA${i} -llapack -lblas -lm 
./TDOA${i} > TDOA${i}.dat

gmt psxy -R0/1500/0/2000 -JX20/10 -Bnesw -Bx100 -By200 -T -K > TDOA${i}.ps
cat TDOA${i}.dat|awk '{if($1=="status:")print $2,$3}'|\
    gmt psxy -R -J -B -Sc0.02c -Gyellow -K -O >> TDOA${i}.ps;
cat TDOA${i}.dat|awk '{if($1=="status:")print $5,$6}'|\
    gmt psxy -R -J -B -Sc0.02c -Gred    -K -O >> TDOA${i}.ps;
gmt psxy -R -J -B -T -O >> TDOA${i}.ps;
gmt psconvert TDOA${i}.ps -A -Tg -E600;
rm TDOA${i} TDOA${i}.ps
done