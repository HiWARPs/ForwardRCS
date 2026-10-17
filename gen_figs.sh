#!/bin/sh
find out -maxdepth 1 -type f | sed s1out/11 | grep -v regge | grep -v _small | grep -v _hw | parallel ./plot.py out/{1} -o fig/{1}.svg -s 9x6
./plot.py out/246_hw -o fig/246_hw.svg -s 5x5
./plot.py out/f_chipt_hw -o fig/f_chipt_hw.svg -s 5x5
./plot.py out/nosr_hw -o fig/nosr_hw.svg -s 5x5
./plot.py out/cs_tot_small -o fig/cs_tot_small.svg -s 5x2.5
./plot.py out/cs_pol_small -o fig/cs_pol_small.svg -s 5x2.5
./plot.py out/re_f_small -o fig/re_f_small.svg -s 3.5x3.5
./plot.py out/f_small -o fig/f_small.svg -s 3x3
./plot.py out/ds_small -o fig/ds_small.svg -s 3.5x2.5
./plot.py out/s2z_small -o fig/s2z_small.svg -s 4x2
