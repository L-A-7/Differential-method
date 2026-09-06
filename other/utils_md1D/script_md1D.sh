#!/bin/sh

# Appel de md1D
#./md1D


#fresults=`./lire_string ./md1D_conf.txt fichier_results`
fresults=$1

ch=/home/lau/Programmes/md1D/utils

# Création de fichiers temporaires
$ch/lire_plot $fresults theta_eff_R eff_R >./eff_R_theta.tmp
$ch/lire_plot $fresults theta_eff_T eff_T >./eff_T_theta.tmp
$ch/lire_plot $fresults profil >./profil.tmp

# Affichage
gnuplot '/home/lau/Programmes/md1D/utils/plotscript'

# Effacement des fichiers .tmp
rm ./eff_R_theta.tmp
rm ./eff_T_theta.tmp
rm ./profil.tmp
