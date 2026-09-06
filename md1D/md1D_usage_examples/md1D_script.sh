#! /bin/sh

ch="/home/lau/Programmes/md1D"
L=16000
h=200
N=100
N_profil=1024
POLA="TE"
# Calcul

			# Création du profil
$ch/utils/profilGen CARRE01 -N_profil $N_profil > ./carre.tmp
$ch/utils/profilGen ALEAT01 -N_profil $N_profil > ./aleat.tmp
$ch/utils/profilGen ADD -f1 ./carre.tmp -f2 ./aleat.tmp -N_profil $N_profil  > ./prof.tmp


echo 'plot "prof.tmp" using 0:1 with linespoints' > ./plot.tmp 
echo 'pause mouse' >> ./plot.tmp
gnuplot ./plot.tmp

prname=ReseauRugueux
$ch/md1D -type_calcul STD -L $L -N $N -h $h -fichier_profil ./prof.tmp \
	-nom_profil $prname -param ./param_ResRug.txt

#for name in {GG5,AY4,REG1}
#do
#	i=1
#	f="$basedir/results/XpC/$name.txt"
#	while [ "$i" -le 10 ]
#	do
#		N1=$((($i-1)*N+1))
#		N2=$(($i*N))
#		echo "type_profil = H_X" > ./tmp.txt
#		echo "N_profil = $N" >> ./tmp.txt
#		echo "profil =" >> ./tmp.txt
#		prname=$name\_TM\_N80_0$i
#		cat $f | $basedir/utils/lire_N1N2 $N1 $N2 >> ./tmp.txt
#		$basedir/md1D -N 80 -fichier_profil ./tmp.txt -nom_profil $prname -param $basedir/md1D_param.txt
#
#
#		i=$(($i+1))
#	done
#done

#for f in $basedir/TESTS/GG* ; do
#	$basedir/md1D -fichier_profil $f -nom_profil $i -param $basedir/md1D_param.txt
#	#g= echo $f | sed -e 's/\/home\/lau\/Programmes\/md1D\/TESTS\///g'
#done

#prname="ech800_400_TE.txt"
#paramfile="param_SPAMM2.txt"
#fichierprofil="$basedir/TESTS/spamm2_L100_CD50_N1024.txt"
#		$basedir/md1D -fichier_profil $fichierprofil -nom_profil $prname -param $paramfile
