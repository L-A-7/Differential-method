#! /bin/sh
# GaelleDASSAULT

dir=$HOME"/Programmes/md1D"
fichier_param="param_RugVol.txt"
N=25
NS=10
L=7000
h=1000


for dh in 1 # {1,0.5,0.2}
do
  for f in vol_binaire1024*  
  do

    prname=results\_$f
    $dir/md1D -N $N -NS $NS -L $L -h $h -fichier_profil $f \
       -delta_h $dh -nom_profil $prname -param $fichier_param
  done
done

