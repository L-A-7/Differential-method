#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

double EcartQ(double period, double CD, double hPoly03, double hThermOx);

int main(){


period=310.0;
CD=130.0;
hPoly03=190.0;
hThermOx=2.0;

N=12
N_profil=1024
angle_i=66.05
delta_h = 2.0

EcartQ(period,CD,hPoly03,hThermox);

	return 0;
}


double EcartQ(double period, double CD, double hPoly03, double hThermOx)
{
	char path ="/home/lau/Programmes/RCWA";
	char options[1000], command[1000];
	int N_profil = 2048;
	double CD, h, L, hPoly03, hThermOx, ecart;

ecart=0.0;

	h = hPoly03 + hThermOx;

	/* Profile generation */
	sprintf(options," CARRE03 -N_profil %d -h1 %d -h2 %d -h3 0 -L1 %d",N_profil,hPoly03/h,hThermOx/h,CD/L);
	sprintf(command,"%f/utils/profilGen %s",path,options);
	system(command);

	/* Differential Method Calculation */
	sprintf(command,"%f/md2D -param param_optim.txt",path);
	system(command);


	/* Reading result from Simulation */

	/* Reading experimental results */
	
	/* Quadratic square difference calculation */

	return ecart;
}
