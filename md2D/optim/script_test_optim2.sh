#!/bin/sh


SCRIPTNAME="./script_fun.sh"
RESFILE="result_file.txt"
PARAMFILE="optim_param.txt"
N_Param=2
N_datas=114


./optim2 -fs $SCRIPTNAME -resfile $RESFILE -paramfile $PARAMFILE -N_param $N_Param -N_datas $N_datas
