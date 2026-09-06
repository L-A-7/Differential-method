#! /bin/sh

dir="./"
dir2=$dir"/utils"
result="result"

param="param_surf_rug.txt"
ref="./surf_rug_ref.txt"

#param="param_compareTayeb.txt"
#ref="./compareTayeb_ref.txt"

file_tmp1="./tmp6549687.tmp"
file_tmp2="./tmp6843219.tmp"

# Calcul 
$dir/md1D -param $param -nom_profil $result

# Creation d'un fichier avec colonnes eff_T et eff_R de chaque fichier resultat 
cat $result.txt | $dir2/lire_tab effR_s > $file_tmp1
cat $ref        | $dir2/lire_tab effR_s | $dir2/add_col $file_tmp1 > $file_tmp2
#cat $result.txt | $dir2/lire_tab eff_R | $dir2/add_col $file_tmp2 > $file_tmp1
#cat $ref        | $dir2/lire_tab eff_R | $dir2/add_col $file_tmp1 > $file_tmp2

# Somme des differences entre les valeurs
echo "Différence par rapport au fichier de référence, en % :"
cat $file_tmp2 | $dir2/diffsum


rm $file_tmp1
rm $file_tmp2

