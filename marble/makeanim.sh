#!/bin/bash

HEIGHT=512
WIDTH=$((${HEIGHT}*16/9))

rm -f frame-*.png

cnt=0
for ppm in `ls -1 marble.chkp-*.ppm`
do
  basename=`echo "$ppm" | cut -d'.' -f1`
  filename=`printf 'frame-%06d.png' ${cnt}`
  printf 'Processing %s ...' ${ppm}
  convert -background black $ppm -gravity center -extent ${WIDTH}x${HEIGHT} ${filename}
  printf 'Done\n'
  cnt=$((cnt+1))
done

ffmpeg  -framerate 1 -i frame-%06d.png -c:v libx264 -y -pix_fmt yuv420p  -vf scale=910:512 -r 30 marble.mp4

