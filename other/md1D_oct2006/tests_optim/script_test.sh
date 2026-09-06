#!/bin/sh

path=$HOME"/Programmes/md1D_oct2006"

period=310.0
CD=126.2
hPoly03=169.3
hThermOx=2.05
N_profil=2048

echo "n1 = 1.5826 + i0.0" > ./tmp.txt
echo "n2 = 1.5826 + i0.0" >> ./tmp.txt

$path/utils/profilGen CARRE03 -N_profil $N_profil -h1 0.988036 -h2 0.011964 -h3 0 -L1 0.407097 >> tmp.txt


$path/md1D -param param_opti.txt -fichier_profil tmp.txt -h 171.35
