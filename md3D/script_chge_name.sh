#! /bin/sh

for f in md* std_in* Makefile
do
	sed 's/md3D_chrono/md_chrono/g' $f > tmp.tmp
	cp tmp.tmp $f
done
