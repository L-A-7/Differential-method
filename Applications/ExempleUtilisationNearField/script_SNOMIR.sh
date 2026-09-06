#! /bin/sh

#PBS -N SNOMIR_TE
#PBS -r n
#PBS -q long
#PBS -l select=1:ncpus=12:mpiprocs=12
#PBS -l walltime=240:00:00

dir=$HOME"/Programs/new_DM"
working_dir=$HOME"/Results/NearField_CopperInSilicon_IR"
param="param_snomIR_edges.txt"
profile_file_without_n=$working_dir"/profile_without_n.txt"
profile_file=$working_dir"/profile.txt"
utils=$HOME"/Programs/utils"
file_ReEx="./ReEx.txt"
file_ReEy="./ReEy.txt"
file_ReEz="./ReEz.txt"
file_ReHpx="./ReHpx.txt"
file_ReHpy="./ReHpy.txt"
file_ReHpz="./ReHpz.txt"
file_ImEx="./ImEx.txt"
file_ImEy="./ImEy.txt"
file_ImEz="./ImEz.txt"
file_ImHpx="./ImHpx.txt"
file_ImHpy="./ImHpy.txt"
file_ImHpz="./ImHpz.txt"
#R_FILE="near_field_results.txt"

export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/opt/intelruntime/11.0.069/

L=2
N_profil=$(echo '512*1' |bc -l)
lambda=10.6 # in microns
h_Cu=0.15
d_Cu=0.25
hmin=-0.1
hmax=0.5
theta_i=45.0
N=100
NS=256
pola=TM
n_Cu='9.128 +i66.433'
#n_Si='3.39 +i0.0'
cd $working_dir

h=$(echo $hmax'-('$hmin')' | bc -l)
#echo "h="$h

for dummy in 50
do
			# Creating raw profile (profile without n)
			d_relatif=$(echo $d_Cu'/'$L | bc -l)
			$dir/../utils/profilGen CARRE02_B -L1 $d_relatif -h1 $h_Cu -h2 0 -N_profil $N_profil > $profile_file_without_n
			# Adding the index value to the profile file
			lambda_Angstrom=$(echo $lambda'*10000' | bc -l)
#			n_Si=$($dir/../utils/refractive_index $HOME/Measurements/Dispersion/Si.dat $lambda_Angstrom)
#			n_Cu=$($dir/../utils/refractive_index $HOME/Measurements/Dispersion/Cu.dat $lambda_Angstrom)

			echo 'n1 = ' $n_Cu > $profile_file
			echo 'h_min = ' $hmin >> $profile_file
			echo 'h_max = ' $hmax >> $profile_file
			cat  $profile_file_without_n >> $profile_file
#			echo 'n_sub = ' $n_Si >> $param

			name="snomIR_edges_L"$L"_h"$h"_dCu"$d_Cu"_hCu"$h_Cu"_"$pola"_N"$N"_NS"$NS"_"$lambda"microns"
			R_FILE="local_field_snomIR_edges_L"$L"_h"$h"_dCu"$d_Cu"_hCu"$h_Cu"_"$pola"_N"$N"_NS"$NS"_"$lambda"microns.txt"
#			echo $name
			$dir/md2D -param $param -N $N -NS $NS -pola $pola -profile_name $name -profile_file $profile_file -h $h -theta_i $theta_i -lambda $lambda -L $L -verbosity 0 > $R_FILE
#			$dir/md2D -param $param -N $N -NS $NS -pola $pola -profile_name $name -profile_file $profile_file -h $h -theta_i $theta_i -lambda $lambda -L $L -calcul_type 'STD'

# DEBUG ######
param_tmp='param_tmp.txt'
echo 'N = ' $N > $param_tmp
echo 'NS = ' $NS >> $param_tmp
echo 'pola = ' $pola >> $param_tmp
echo 'h = ' $h >> $param_tmp
echo 'L = ' $L >> $param_tmp
echo 'theta_i = ' $theta_i >> $param_tmp
echo 'lambda = ' $lambda >> $param_tmp
echo 'profile_name = ' $name >> $param_tmp
echo 'profile_file = ' $profile_file >> $param_tmp
echo 'verbosity = 0' >> $param_tmp
cat $param >> $param_tmp

##############


done
cat $R_FILE | $utils/lire_tab2 "sumT11"

cat $R_FILE | $utils/lire_tab2 "ReEx" >  $file_ReEx
cat $R_FILE | $utils/lire_tab2 "ReEy" >  $file_ReEy
cat $R_FILE | $utils/lire_tab2 "ReEz" >  $file_ReEz
cat $R_FILE | $utils/lire_tab2 "ReHpx" > $file_ReHpx
cat $R_FILE | $utils/lire_tab2 "ReHpy" > $file_ReHpy
cat $R_FILE | $utils/lire_tab2 "ReHpz" > $file_ReHpz
cat $R_FILE | $utils/lire_tab2 "ImEx" >  $file_ImEx
cat $R_FILE | $utils/lire_tab2 "ImEy" >  $file_ImEy
cat $R_FILE | $utils/lire_tab2 "ImEz" >  $file_ImEz
cat $R_FILE | $utils/lire_tab2 "ImHpx" > $file_ImHpx
cat $R_FILE | $utils/lire_tab2 "ImHpy" > $file_ImHpy
cat $R_FILE | $utils/lire_tab2 "ImHpz" > $file_ImHpz


