----------------------- md1D README.txt -----------------------

md1D is a program calculating the diffraction of an electromagnetic wave by a structure.
md1D stands for "Methode Differentielle 1D" (One Dimension Differential Method)

INSTALLATION

installation needs the libraries fftw and GSL, then simply type 'make' in the directory containing the sources

USAGE

For basic usage, type 'dir/md1D' in the console, 'dir' being the directory where is the executable 'md1D'
By default the parameters are set in the config file 'md1D_param.txt', it's generaly more convenient
to set them in an other file (one file for one type of calculus) and to call md1D with
'md1D -param param_file_name'.
The profile must be defined in a separate file, its header containes a few informations on the profile
and the rest of the file contains the profile definition datas.

For a more advanced usage, see the scripts, parameters files and profile files given as examples

