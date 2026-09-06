#! /bin/sh
# Rugo surface

dir=$HOME"/Programmes/md1D"
N=400
NS=20
Nx=1024
angle_i=50
delta_h=1.0


for f in rug1D_R0100*
do
  L=$(echo $f | sed 's/....................._L//g' | sed 's/_Nx.....txt//g')
  Nx=$(echo $f | sed 's/.....................................//g' | sed 's/.txt//g')

echo "L="$L" Nx="$Nx

  # cat $f | sed 's/  /\n/g' > tmp.txt
  # mv tmp.txt $f

  echo "type_profil=H_X" > ./tmp.txt
  echo "N_x = $Nx" >> ./tmp.txt
  echo "profil =" >> ./tmp.txt
  cat $f >> ./tmp.txt

	prname=results$f
	
	$dir/md1D -L $L -N $N -NS $NS -fichier_profil ./tmp.txt \
	    -delta_h $delta_h -angle_i $angle_i -nom_profil $prname -param ./param_verresable_i50.txt

done

