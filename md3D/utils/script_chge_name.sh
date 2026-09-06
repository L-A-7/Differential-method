#! /bin/sh

for f in ./*.c 
do
	sed 's/md2D/md3D/g' $f > tmp.tmp
	cp tmp.tmp $f
done
