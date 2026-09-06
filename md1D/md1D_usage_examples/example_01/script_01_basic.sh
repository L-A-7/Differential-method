#! /bin/sh
# Script showing a simple usage of md1D
# most parameters are configured in param_01_basic.txt

# Path to the program
path=$HOME"/Programmes/md1D"

# Calling md1D with 2 options
$path/md1D -param ./param_01_basic.txt

# -param param_file : indicate a particular param_file for md1D, 
#                     by default, md1D looks in param_md1D.txt if the file exists
