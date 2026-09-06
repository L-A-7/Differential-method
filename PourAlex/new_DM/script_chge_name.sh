#! /bin/sh

for f in md* std_in* Makefile
do
	sed 's/fichier_profil/profile_file/g' $f > tmp.tmp
	cp tmp.tmp $f
done
