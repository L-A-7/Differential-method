#include "optimiz.h"

int main(){

/*	double period=310.0;*/
	char command[5000];
	double hPoly03=170.0;
	double hThermOx=2.0;
	double period = 310.0;
	double param[10];
	int Nb_lambda;
	double dNb_lambda, *lambda_exp, *tanPsicosDelta_exp, *tanPsicosDelta_simu, *tanPsi_exp, *cosDelta_exp;
	double CDmin, CDmax, CDpas, h1min, h1max, h1pas, h1, CD, EqQuadra, EqQuadra_tanPsi, EqQuadra_cosDelta;
	int NCD, Nh1,nl;

	/* Reading experimental results */
	lire_tab("./measured_data.txt", "Nb_lambda", &dNb_lambda, 1);
	Nb_lambda = ROUND(dNb_lambda);
	lambda_exp = (double *) malloc(sizeof(double)*Nb_lambda);
	tanPsicosDelta_exp  = (double *) malloc(sizeof(double)*2*Nb_lambda);
	tanPsicosDelta_simu = (double *) malloc(sizeof(double)*2*Nb_lambda);
	cosDelta_exp = &tanPsicosDelta_exp[Nb_lambda];
	tanPsi_exp   = &tanPsicosDelta_exp[0];
	lire_tab("./measured_data.txt", "var_lambda",  lambda_exp, Nb_lambda);
	lire_tab("./measured_data.txt", "var_tan_Psi",  tanPsi_exp, Nb_lambda);
	lire_tab("./measured_data.txt", "var_cos_delta",  cosDelta_exp, Nb_lambda);

	CDmin = 114;
	CDmax = 116;
	NCD = 5;
	h1min = 178;
	h1max = 182;
	Nh1 = 9;
	
	CDpas=(CDmax-CDmin)/(NCD-1);
	h1pas=(h1max-h1min)/(Nh1-1);

	sprintf(command,"echo \"CD hPoly EcartQ_tanPsi EcartQ_cosDelta\" > ./carto.txt"); system(command);
	
	for (CD=CDmin;CD<=CDmax;CD+=CDpas){
		for (h1=h1min;h1<=h1max;h1+=h1pas){
			hPoly03 = h1;
			param[0] = CD;
			param[1] = hPoly03;
			param[2] = hThermOx;
			param[3] = period;
		
			tanPsi_cosDelta_lambda(tanPsicosDelta_simu, lambda_exp, Nb_lambda, param);
			EqQuadra_tanPsi = 0;
			EqQuadra_cosDelta = 0;
			for (nl=0;nl<=Nb_lambda;nl++){
				EqQuadra_tanPsi   += (tanPsicosDelta_simu[nl] - tanPsi_exp[nl])*(tanPsicosDelta_simu[nl] - tanPsi_exp[nl])/Nb_lambda;
				EqQuadra_cosDelta += (tanPsicosDelta_simu[nl+Nb_lambda] - cosDelta_exp[nl])*(tanPsicosDelta_simu[nl+Nb_lambda] - cosDelta_exp[nl])/Nb_lambda;
			}
		sprintf(command,"echo \"%f %f %f %f\" >> ./carto.txt",CD,hPoly03,EqQuadra_tanPsi, EqQuadra_cosDelta); system(command);
		}
	}

	free(tanPsicosDelta_simu);
	free(tanPsicosDelta_exp);
	free(lambda_exp);

	return 0;
}



double fonction_simple(double x)
{

	return cos(x*3.141592/180.0);
}

void func(double *parameters, double *tanPsicosDelta, int m, int n, void *adata)
{
	int i, Nb_lambda, param_out_of_range;
	double period, CD, hPoly03, hThermOx, *var_lambda;
	double param[10];

	CD = parameters[0];
	hPoly03 = parameters[1];
	hThermOx = parameters[2];
	period = 310.0;

	param[0] = CD;
	param[1] = hPoly03;
	param[2] = hThermOx;
	param[3] = period;

	Nb_lambda = n/2;
	var_lambda = (double *) adata;

	/* verification des limites */
	param_out_of_range = 0;
	if (CD<=0 || CD >= period || hPoly03 <= 0 || hThermOx <= 0 || hThermOx >= 10){
		param_out_of_range = 1;
	}
	
	if (param_out_of_range){
		for(i=0;i<=2*Nb_lambda-1;i++){
			tanPsicosDelta[i] = 1.0e10; /* results put to ~infinity, i.e. bad point */
		}
	return;
	}

	tanPsi_cosDelta_lambda(tanPsicosDelta, var_lambda, Nb_lambda, param);

}


int tanPsi_cosDelta_lambda(double *tanPsicosDelta, double *var_lambda, int Nb_lambda, double *param)
{
	FILE *fp;
	char *path ="/home/lau/Programmes/RCWA";
	char options[5000], command[5000], results_file[5000];
	int N_profil = 2048;
	int nl, i;
	double *var_tan_Psi, *var_cos_delta;
	complex index[10];
	
	double CD = param[0];
	double hPoly03 = param[1];
	double hThermOx = param[2];
	double period = param[3];
	char material[10][1000];
	int N = 8;
	int NS;
	int N_layers = 2;
	char tab_NS_string[5000];
	int ns_hPoly, NS_hPoly;
	double *tab_NS;
	
	tab_NS = (double *) malloc(sizeof(double)*1000);
	
	var_tan_Psi   = tanPsicosDelta;
	var_cos_delta = &tanPsicosDelta[Nb_lambda];

	sprintf(material[0],"AIR");
	sprintf(material[1],"POLY03");
	sprintf(material[2],"OXIDE_THERM");
	sprintf(material[3],"SI_CRISTAL");

	/* Boucle sur lambda */
	for (nl=0;nl<=Nb_lambda-1;nl++){

		/* lambda determination */
		printf("lambda = %f\n",var_lambda[nl]);
		
		/* Refractive indices determination */
		for (i=0; i<=N_layers+1; i++){
			index[i] = refractive_index(material[i],var_lambda[nl],"Lookup");
		}

	/* tab_NS making */ /* precise le decoupage en tranches matrice S */
	tab_NS[0]=0;
	tab_NS[1]=hThermOx;
	NS_hPoly = CEIL(10*hPoly03/var_lambda[nl]);
	NS=1;
	sprintf(tab_NS_string,"tab_NS = %f %f",tab_NS[0],tab_NS[1]);
	for (ns_hPoly=1;ns_hPoly<=NS_hPoly;ns_hPoly++){
		tab_NS[ns_hPoly+1]=tab_NS[1]+ns_hPoly*(hPoly03/NS_hPoly);
		sprintf(tab_NS_string,"%s %f",tab_NS_string,tab_NS[ns_hPoly+1]);
		NS++;
	}
	/*printf("%s\n",tab_NS_string);*/
	sprintf(command,"echo \"%s\" > ./tab_NS.txt",tab_NS_string); system(command);

		/* Profile generation */
		sprintf(command,"echo \"n1 = %f + i%f\" > ./tmp.txt",creal(index[1]),cimag(index[1])); system(command);
		sprintf(command,"echo \"n2 = %f + i%f\" >> ./tmp.txt",creal(index[2]),cimag(index[2])); system(command);
		sprintf(options," CARRE03 -N_profil %d -h1 %f -h2 %f -h3 0 -L1 %f",N_profil,hPoly03,hThermOx,CD/period);
		sprintf(command,"%s/utils/profilGen %s >> tmp.txt",path,options);
		system(command);

		/* param_file making */
		sprintf(command,"echo \"n_super = %f + i%f\" > ./param_opti.txt",creal(index[0]),cimag(index[0])); system(command);
		sprintf(command,"echo \"n_sub = %f + i%f\" >> ./param_opti.txt",creal(index[3]),cimag(index[3])); system(command);
		sprintf(command,"echo \"h = AUTO\" >> ./param_opti.txt"); system(command);
		sprintf(command,"echo \"i_field_mode  = PLANE_WAVE\" >> ./param_opti.txt"); system(command);
		sprintf(command,"echo \"calcul_type   = ELLIPSO\" >> ./param_opti.txt"); system(command);
		sprintf(command,"echo \"calcul_method = RCWA\" >> ./param_opti.txt"); system(command);
		sprintf(command,"echo \"delta_h = 0.1\" >> ./param_opti.txt"); system(command);


		/* Calling the Calculation Method */
		sprintf(command,"%s/md2D -param param_opti.txt -nom_profil one_lambda_results -fichier_profil tmp.txt",path);
		sprintf(command,"%s -L %f -lambda %f -N %d -NS %d -theta_i 66.05 -psi 0 -phi_i 0 -verbosity 1",command,period,var_lambda[nl],N,NS);
		sprintf(command,"%s -tab_NS_ENABLED 1 -tab_NS_filename tab_NS.txt",command);
		system(command);

		/* Reading tan_Psi & cos_Delta */
		lire_tab("./one_lambda_results.txt", "spec_tan_Psi", &var_tan_Psi[nl], 1);
		lire_tab("./one_lambda_results.txt", "spec_cos_delta", &var_cos_delta[nl], 1);
	}

	/* writing results to a file */
	sprintf(results_file,"results_CD%f_hPoly%f.txt",CD,hPoly03);
	if (!(fp = fopen(results_file,"w+"))){
		fprintf(stderr, "%s ligne %d : ERREUR, impossible d'ouvrir %s\n",__FILE__, __LINE__,"results_var_lambda.txt");
		exit(EXIT_FAILURE);
	}
	fprintf(fp, "var_lambda    = ");	ecrire_dble_tab(fp, var_lambda, Nb_lambda, " ", 250, "\n");
	fprintf(fp, "\nvar_tan_Psi   = ");	ecrire_dble_tab(fp, var_tan_Psi, Nb_lambda, " ", 250, "\n");
	fprintf(fp, "\nvar_cos_delta = ");	ecrire_dble_tab(fp, var_cos_delta, Nb_lambda, " ", 250, "\n");
	fclose(fp);
	free(tab_NS);
/*system("scite ./results_var_lambda.txt &");*/
	return 0;
}



/*---------------------------------------------------------------------------------------------*/
/*!	\fn		int lire_tab(char *nom_fichier, const char *label, double *tab, int N)
 *
 *	\brief	Lit N valeurs de format double dans un fichier et les stocke dans un tableau   \n 
 *			Les valeurs doivent être séparées par un ou plusieurs espaces, tabulations     \n
 *			ou sauts de lignes et précédées d'un label éventuellement suivi d'un signe '='.\n
 *          ex. : (...) tab1 = 3.4  4.5e-3  +46  -7.6e+2 ...                               \n
 *          Remarque : Pour lire un tableau sans label, donner "" comme label.             \n
 *
 *	\return	0 si succès, 1 si le nombre d'éléments lus diffère de N ou si le label n'a pas été trouvé.
 *
 *	\todo	RENDRE PLUS ROBUSTE : PAS DE BUFFER OVERFLOW AU CAS OU IL Y A PLUS DE N LIGNES
 *
 */
/*---------------------------------------------------------------------------------------------*/
int lire_tab(const char *nom_fichier, const char *label, double *tab, int N)
{
	FILE *fp;
	char *pos, *pos2, *endptr, line[SIZE_LINE_BUFFER];
	int cpt=0, line_cpt=0;
	double tmp;
	
	if (!(fp = fopen(nom_fichier,"r"))){
		fprintf(stderr, "%s ligne %d : ERREUR, impossible d'ouvrir %s\n",__FILE__, __LINE__,nom_fichier);
		return 1;
	}

	/* Recherche du label */
	while(!feof(fp)){ 
		if (lire_ligne(fp,line) != 0) goto LECTURE_FINIE;
		line_cpt++;
		skip_comment(line);
		if ((pos = label_search(line,label)) != NULL){ /* on cherche le label */	
			pos += strlen(label); 
			if((pos2 = strchr(pos,'=')) != NULL) *pos2 = ' '; /* on remplace '=' par ' ' */
			goto LABEL_TROUVE;
		}
	}

	fprintf(stderr,"%s ligne %d : ERREUR, le label '%s' n'a pas été trouvé dans %s \n", __FILE__, __LINE__, label, nom_fichier);
	fclose(fp);
	return 1;

	/* Lecture des valeurs */
	while(!feof(fp)){ /* Tant qu'on est pas à la fin du fichier */
		if (lire_ligne(fp,line) != 0) {goto LECTURE_FINIE;}
		line_cpt++;
		skip_comment(line); 
		pos = line;
		
	LABEL_TROUVE :
		while(isspace(*pos)) pos++; /* on élimine les espaces */
		while(pos < line+strlen(line)) { /* Tant qu'on est pas à la fin de la ligne */
			tmp = strtod(pos, &endptr); /* on lit le 'double' */
			if (pos == endptr) {goto LECTURE_FINIE;} /* conversion ratée */
			if (cpt <= N+1) tab[cpt] = tmp;
			cpt ++;
			pos = endptr;
			while(isspace(*pos)) pos++; /* on élimine les espaces */
		}
	}

	LECTURE_FINIE:
	if (cpt != N) {
		fprintf(stderr, "%s ligne %d : ERREUR, %s contient %d valeurs au lieu de %d dans %s (ligne %d)\n", __FILE__, __LINE__, label, cpt, N, nom_fichier, line_cpt);
		fclose(fp);
		return 1;
	}
	fclose(fp);
	return 0;
}

/*---------------------------------------------------------------------------------------------*/
/*! \fn		void lire_ligne(FILE *fp, char *line)
 *
 *  \brief	Lit une ligne dans un fichier et la stocke dans une chaine de charactères 
 */
/*---------------------------------------------------------------------------------------------*/
int lire_ligne(FILE *fp, char *line)
{
	if (fgets(line, SIZE_LINE_BUFFER, fp) == NULL) return 1;
	if (strlen(line) == SIZE_LINE_BUFFER-1) {
		fprintf(stderr, "%s ligne %d : ERREUR, taille de buffer insuffisante,impossible de lire plus de "
						"%d caractères par ligne.\n",__FILE__, __LINE__,SIZE_LINE_BUFFER-1);
		exit(EXIT_FAILURE);
	}
	return 0;
}


/*---------------------------------------------------------------------------------------------*/
/*!	\fn		char *label_search(char *str,const char *label)
 *
 *	\brief	Cherche un label dans une chaine de caractères, le label doit être isolé, c.a.d, \n
 *          en début de ligne ou précédé d'un espace au sens de isspace() et suivi d'un espace \n
 *          ou d'un signe '='
 *
 *	\return	La position de la 1ere occurence du label dans la chaine ou NULL si le label n'a pas été trouvé
 */
/*---------------------------------------------------------------------------------------------*/
char *label_search(char *str,const char *label)
{
	char *pos;

	/* Cas du label vide */
	if (strlen(label) == 0) return str;

	/* Recherche du label */
	while((pos=strstr(str,label)) != NULL) {
		/* Vérification que le label est en début de ligne ou précédé par un espace */
		if (pos != str && !isspace(*(pos-1))) {
			str = pos + strlen(label);
			continue;
		}
		/* Vérification que le label est suivi par un espace, saut de ligne ou signe '=' */
		if (strlen(pos) > strlen(label)) {
			if (!isspace(*(pos+strlen(label))) && *(pos+strlen(label)) != '=') {
				str = pos + strlen(label);
				continue;
			}
		}
		break;
	}
	return pos;
}


/*---------------------------------------------------------------------------------------------*/
/*! \fn		void skip_comment(char *str_in_out)
 *
 *  \brief	Elimine tout ce qui se trouve après un commentaire '#' dans str_in_out
 */
/*---------------------------------------------------------------------------------------------*/
void skip_comment(char *str_in_out){

	char *pos;
	/* Cherche CHAR_COMMENT et le remplace par le charactère nul '\0' */
	if((pos = strchr(str_in_out,CHAR_COMMENT)) != NULL) {
		*pos = '\0';
	}
 }

/*---------------------------------------------------------------------------------------------*/
/*!	\fn		int	ecrire_dble_tab(FILE *fp, double *tab, int N, char *separateur1, int Nmax1, char *separateur2)
 *
 *	\brief	Ecrit les valeurs d'un tableau séparées par les 'séparateurs1' (par ex " "), plus par les \n
 *          'séparateurs2' (par ex "\n") une fois tous les Nmax1 éléments.
 */
/*---------------------------------------------------------------------------------------------*/
int ecrire_dble_tab(FILE *fp, double *tab, int N, char *separateur1, int Nmax1, char *separateur2)
{
	int i,k=0;

	while((k+1)*Nmax1 < N) {
		for (i=k*Nmax1; i<=MIN(N-1,(k+1)*Nmax1-1); i++){
			fprintf(fp,"% 1.6e%s",tab[i],separateur1);
		}
		k++;
		fprintf(fp,"%s",separateur2);
	}
	/* Derniere ligne, traitée à part car pas de séparateur2 à la fin ! */
	for (i=k*Nmax1; i<=MIN(N-1,(k+1)*Nmax1-1); i++){
		fprintf(fp,"% 1.6e%s",tab[i],separateur1);
	}
	
	return 0;
}

/*-------------------------------------------------------------------------------------*/
/*!	\fn		int SaveDbleTab2file (char *filename, double *tab, int N, char *separateur)
 *
 *	\brief
 */
/*-------------------------------------------------------------------------------------*/
int SaveDbleTab2file (double *tab, int N, char *filename, char *separateur)
{
	int i;
	FILE *fp;
	
	/* Affichage à l'écran si filename = "stdout" */										
	if(!strcmp(filename,"stdout")) fp = stdout;
	else fp = fopen(filename, "w");
	
	for (i=0;i<=N-1;i++){
		fprintf(fp,"%1.6e%s",tab[i],separateur);
	}

	if(strcmp(filename,"stdout")) fclose(fp);

	return 0;
}


/*---------------------------------------------------------------------------------------------*/
/*!	\fn		complex md2D_index(char* name,double lambda, char* method)
 *
 *	\brief
 */
/*---------------------------------------------------------------------------------------------*/
complex refractive_index(char* name,double lambda, char* method)
{
	char filename[SIZE_STR_BUFFER];
	double lambda_angstrom = 10*lambda;
	double cauchy_n[5],cauchy_k[5],n_power[6],k_power[6],index_n,index_k;
	complex index;
	int i,Nb_cauchy = 5;
	
	/* File name */
	snprintf(filename, SIZE_STR_BUFFER*sizeof(char), "%s_index.txt", name);
		
	/* Cauchy method */
	if (!strcmp(method,"Cauchy")){
		lire_tab(filename, "CAUCHY_N", cauchy_n, Nb_cauchy);
		lire_tab(filename, "CAUCHY_K", cauchy_k, Nb_cauchy);
		lire_tab(filename, "N_POWERS", n_power, Nb_cauchy+1);
		lire_tab(filename, "K_POWERS", k_power, Nb_cauchy+1);
/*SaveDbleTab2file (cauchy_n,  Nb_cauchy,"stdout", " ");printf("\n");
SaveDbleTab2file (n_power, Nb_cauchy,"stdout", " ");printf("\n");
SaveDbleTab2file (cauchy_k,  Nb_cauchy,"stdout", " ");printf("\n");
SaveDbleTab2file (k_power, Nb_cauchy,"stdout", " ");printf("\n");*/

		index_n = 0;
		index_k = 0;
		for (i=0;i<=Nb_cauchy-1;i++){
			index_n += cauchy_n[i]*pow(lambda_angstrom,n_power[i]);
			index_k += cauchy_k[i]*pow(lambda_angstrom,k_power[i]);
		}
/*printf("Milieu: %s, lambda = %f, n = %f, k = %f\n",name,lambda,index_n,index_k);
*/				
		index = index_n + I*index_k; 
		return index;
		
	}
	/* Lookup-Table method */
	if (!strcmp(method,"Lookup")){
		double *tab_n, *tab_k, *tab_lambda, npoints, *table_tmp;
		lire_tab(filename, "NPOINTS", &npoints, 1);
		table_tmp= (double *) malloc(sizeof(double)*3*((int)npoints));
		tab_lambda= (double *) malloc(sizeof(double)*((int)npoints));
		tab_n= (double *) malloc(sizeof(double)*((int)npoints));
		tab_k= (double *) malloc(sizeof(double)*((int)npoints));
		lire_tab(filename, "TABLE", table_tmp, 3*(int)npoints);
		/* Réarrangement en plusieurs tableaux */
		for (i=0; i<=(int)npoints -1;i++){
			tab_lambda[i] = table_tmp[3*i];
			tab_n[i] = table_tmp[3*i+1];
			tab_k[i] = table_tmp[3*i+2];
		}
		/* Recherche de l'indice de tableau pour lambda */
		int num_min = 0;
		int num_max = (int) npoints-1;
		int numero = num_max>>1;
		while(!(tab_lambda[numero] <= lambda && lambda < tab_lambda[numero+1])){
			numero=num_min+((num_max-num_min)>>1);
			if (lambda < tab_lambda[numero]){
				num_max=numero;
			}
			if (tab_lambda[numero+1] <= lambda){
				num_min=numero;
			}
			if (num_min == num_max){
				fprintf(stderr,"md2D_indice, look-up table incompatible avec lambda = %f pour %s",lambda,name);
				exit(EXIT_FAILURE);
			}
		}
		/* Interpolation linéaire */
		double l1 = tab_lambda[numero];
		double l2 = tab_lambda[numero+1];
		double n1 = tab_n[numero];
		double n2 = tab_n[numero+1];
		double k1 = tab_k[numero];
		double k2 = tab_k[numero+1];
		index_n= (lambda-l1)*(n2-n1)/(l2-l1) + n1;
		index_k= (lambda-l1)*(k2-k1)/(l2-l1) + k1;

		free(tab_lambda);
		free(tab_n);
		free(tab_k);
		free(table_tmp);
		
/*printf("lambda= %f, lambda1= %f, lambda2= %f\n",lambda,tab_lambda[numero],tab_lambda[numero+1]);
printf("lambda= %f, n = %f, k= %f\n",lambda,index_n,index_k);
*/		index = index_n + I*index_k; 

		return index;
	}
	return -1;/* ERREUR, ne devrait pas arriver là ... */
}
