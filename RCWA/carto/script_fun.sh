#!/bin/sh

# Reconstruction d'un profil carre avec couche d'oxide (cf profilGen.c, CARRE03)
# utilisation : script_fun.sh CD hPoly hThermOx

ch=$HOME"/Programmes/RCWA"
util=$ch"/optim"
PROFIL="./profil_tmp.txt"
S_STEP_FILE="imposed_S_step.txt"
RESULT_FILE="result_file.txt"
SIMU_DATA_FILE="simu_data.txt"
EXP_DATA_FILE="measured_data01.txt"
param="param_reconst.txt"
lambda_min=210
lambda_Max=770
Period=310.0
N=15
NS=10
# Reading arguments
CD=$1
hPoly=$2
hThermOx=2.0
echo "CD = "$CD", hPoly = "$hPoly 
>tan_Psi.tmp
>cos_delta.tmp
>lambda.tmp

# checking arg values
out_of_limits=0
if [ $(echo "$CD > $Period"|bc) -eq 1 ]; then
	out_of_limits=1
elif [ $(echo "$CD < 0 "|bc) -eq 1 ]; then
	out_of_limits=1
elif [ $(echo "$hPoly < 0 "|bc) -eq 1 ]; then
	out_of_limits=1
elif [ $(echo "$hThermOx < 0 "|bc) -eq 1 ]; then
	out_of_limits=1
fi
if [ $out_of_limits -gt 0 ]; then # return pseudo infinity if out of limits
	cat ./big_values_114.txt > $RESULT_FILE
	exit
fi

# imposed S steps
echo $hThermOx > $S_STEP_FILE

# loop over lambda
lambda=$lambda_min
Nlambda=0
while [ "$lambda" -le $lambda_Max ]
do
#echo $lambda
	# Refractive index determination
	n0=$($util/refractive_index AIR $lambda)
	n1=$($util/refractive_index POLY03 $lambda)
	n2=$($util/refractive_index OXIDE_THERM $lambda)
	n3=$($util/refractive_index SI_CRISTAL $lambda)

	# Creating the profile
	CD_ratio=$(echo $CD/$Period | bc -l)
	echo "n1 = "$n1 > $PROFIL
	echo "n2 = "$n2 >> $PROFIL
	$ch/utils/profilGen CARRE03 -N_profil 131072 -h1 $hPoly -h2 $hThermOx -h3 0 -L1 $CD_ratio >> $PROFIL
	PROFIL_NAME="carre_CD"$CD"_hPoly"$hPoly"_hThermOx"$hThermOx

	# Refractive index values to param_file
	cat $param | grep --invert-match n_su > param.tmp
	echo "n_super = "$n0 > $param
	echo "n_sub   = "$n3 >> $param
	cat param.tmp >> $param

	 
	# calling the differential method
	$ch/md2D -imposed_S_steps 1 -imposed_S_steps_filename $S_STEP_FILE -N_imposed_S_steps 1 \
				-param $param -L $Period -fichier_profil $PROFIL -nom_profil $PROFIL_NAME -verbosity 0 \
				-lambda $lambda -N $N -NS $NS

	# extracting results to temporary files
	cat $PROFIL_NAME.txt | $ch/utils/lire_tab spec_tan_Psi >> tan_Psi.tmp
	cat $PROFIL_NAME.txt | $ch/utils/lire_tab spec_cos_delta >> cos_delta.tmp
	echo $lambda >> lambda.tmp
	lambda=$(($lambda+10))
	Nlambda=$(($Nlambda+1))
done

# Writting results to result_file
echo "lambda =" > $SIMU_DATA_FILE
cat lambda.tmp | $ch/utils/lire_tab "" >> $SIMU_DATA_FILE
echo "tan_Psi =" >> $SIMU_DATA_FILE
cat tan_Psi.tmp | $ch/utils/lire_tab "" >> $SIMU_DATA_FILE
echo "cos_Delta =" >> $SIMU_DATA_FILE
cat cos_delta.tmp | $ch/utils/lire_tab "" >> $SIMU_DATA_FILE

cp $SIMU_DATA_FILE "res_"$PROFIL_NAME".txt"

# making substraction with measured data and checking that the same lambda are used
$ch/optim/sub_tanpsicosdelta $Nlambda $SIMU_DATA_FILE $EXP_DATA_FILE > $RESULT_FILE

