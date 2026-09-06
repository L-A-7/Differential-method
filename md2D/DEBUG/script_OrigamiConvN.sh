#! /bin/sh
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/opt/intel/compilers_and_libraries_2018.2.199/linux/mkl/lib/intel64_lin/
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/opt/intel/compilers_and_libraries_2018.2.199/linux/compiler/lib/intel64_lin/

exec_dir=$HOME"/Programs/md2D"
utils_dir=$HOME"/Programs/utils"
working_dir=$HOME"/Programs/md2D/DEBUG"

param="param_Origami.txt"
profile_file_without_n=$working_dir"/profile_without_n.txt"
profile_file=$working_dir"/profile.txt"

L=1000
N_profil=$(echo '4096*1' |bc -l) 
theta_i=0
N=50
NS=1000
pola=TM
alpha_deg=54.7
PI=3.14159265358979323846
h1=50
n2="1.38 +i0.0" # index of glue substrate

tan_alpha=$(echo "s($alpha_deg*$PI/180)/c($alpha_deg*$PI/180)"|bc -l)
h2=$(echo "0.5*$tan_alpha*$L"|bc -l) 

h=$(echo "$h1+$h2"| bc -l)


echo $h


cd $working_dir

for N in 10 20 30 40 50 60 70 80 90 100 110 120 130 140 150 #50 40 60 30 20 70
do
	for lambda in 600 #400 800 #
	do
	# Creating raw profile (profile without n)
	$utils_dir/profilGen ORIGAMI -h1 $h1 -h2 $h2 > $profile_file_without_n
	# Adding the index value to the profile file
	lambda_Angstrom=$(echo $lambda'*10' | bc -l)
	n1=$($utils_dir/refractive_index $HOME/Measurements/Dispersion/Au.dat $lambda_Angstrom)
	echo 'n1 = ' $n1
	echo 'n1 = ' $n1 > $profile_file
	echo 'n2 = ' $n2
	echo 'n2 = ' $n2 >> $profile_file
	cat  $profile_file_without_n >> $profile_file

	name="W_test_N"$N"_NS"$NS"_"$lambda"nm_"$pola
	echo $name
	$exec_dir/md2D -param $param -N $N -NS $NS -pola $pola -profile_name $name -profile_file $profile_file -h $h -theta_i $theta_i -lambda $lambda -L $L
done
done



	

