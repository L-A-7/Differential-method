/*---------------------------------------------------------------------------------------------*/
/*!	\file		md1D_pilot.c
 *
 * 	\brief		Pilotage de md1D 
 */
/*---------------------------------------------------------------------------------------------*/


#include "md1D_pilot.h"
/*---------------------------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------------------------*/
int main(int argc, char *argv[]){

	struct Noms_fichiers nomfichier;
	struct Param_struct param;
	struct Efficacites_struct effic;

	param.verbose =1;

	/* Copie des arguments de la ligne de commande */
	int i;
	char **argvcp; 
	argvcp = (char **) malloc(sizeof(char*)*argc);
	argvcp[0] = (char*) malloc(sizeof(char)*SIZE_STR_BUFFER*argc);
	for(i=1;i<=argc-1;i++){
		argvcp[i] = argvcp[i-1] + SIZE_STR_BUFFER;
		strncpy(argvcp[i], argv[i],SIZE_STR_BUFFER);
	}
	param.argc = argc;
	param.argvcp = argvcp;

	/* Initialisation du programme : lecture des données, allocation de mémoire, etc. */
	md1D_init(&param, &effic, &nomfichier);

/*
complex ***T;
int j, k, ntab=32, ncol=31, nlign=31;
T = allocate_CplxMatrix_3(nlign,ncol,ntab);
for(k=0;k<ntab;k++){
	for(i=0;i<ncol;i++){
		for(j=0;j<nlign;j++){
			T[k][i][j] = 0.0*i + 0.0*j + 1.0*k;}}}

SaveMatrix2file (T[31], nlign, ncol, "Re","stdout");printf("\n");

return 0;
*/
		
	/* Choix du type de calcul */
	if (!strcmp(param.calcul_type,"STD")){
		md1D_standard (&param, &effic, &nomfichier);
	}else if (!strcmp(param.calcul_type,"VAR_I")){
		md1D_variation_incidence_B (&param, &effic, &nomfichier);
	}else if (!strcmp(param.calcul_type,"VAR_I_BIS")){
		md1D_variation_incidence_A (&param, &effic, &nomfichier);
	}else if (!strcmp(param.calcul_type,"ELLIPSO")){
		md1D_ellipso (&param, &effic, &nomfichier);
	}else if (!strcmp(param.calcul_type,"VAR_I_ELLIPSO")){
		md1D_var_i_ellipso (&param, &effic, &nomfichier);
	}else if (!strcmp(param.calcul_type,"VAR_LAMBDA_ELLIPSO")){
		md1D_var_lambda_ellipso (&param, &effic, &nomfichier);
	}else if (!strcmp(param.calcul_type,"ALEAT_T_ELLIPSO")){
		md1D_aleat_T_ellipso (&param, &effic, &nomfichier);
	}else if (!strcmp(param.calcul_type,"REPLIC_T_ELLIPSO")){
		md1D_replic_T_ellipso (&param, &effic, &nomfichier);
	}else if (!strcmp(param.calcul_type,"NEAR_FIELD")){
		md1D_near_field (&param, &effic, &nomfichier);
	}else{
		fprintf(stderr, "%s, ligne %d : Erreur, type de calcul inconnu (\"%s\")\n",__FILE__,__LINE__,param.calcul_type);
		exit(EXIT_FAILURE);
	}
		

printf("Somme des efficacités   : %1.10f\n",effic.somm_eff);
printf("1-Somme des efficacités : %e\n",1-effic.somm_eff);


/*************************************** TESTS *********************************************************/
/*complex *FFT_invk2_directe(double z, complex *TF_invk2, struct Param_struct *par);
complex *FFT_k2_directe(double z, complex *TF_k2, struct Param_struct *par);
complex *FFT_k2_stockee(double z, complex *TF_k2, struct Param_struct *par);

int i,j, Nb=300;
double z;
complex *TF_invk2  = (complex *) malloc(sizeof(complex)*(4*param.N+1));
complex *TF_k2     = (complex *) malloc(sizeof(complex)*(4*param.N+1));
TF_invk2 = FFT_invk2_directe(0.5*param.h, TF_invk2, &param);

SaveCplxTab2file (TF_invk2, 4*param.N+1, "Re", "stdout");*/
/*int N = param.N;*/
/*SaveMatrix2file (param.S12, 2*N+1, 2*N+1, "Re", "stdout");
SaveMatrix2file (param.S22, 2*N+1, 2*N+1, "Re", "stdout");
SaveMatrix2file (param.S12, 2*N+1, 2*N+1, "Im", "stdout");
SaveMatrix2file (param.S22, 2*N+1, 2*N+1, "Im", "stdout");*/

/************************************* FIN TESTS *******************************************************/


	/* Libération de la mémoire */
	md1D_free(&param, &effic);
	
	/* Affichage du temps écoulé */
	fprintf(stdout,"Temps écoulé : %f s \n", md1D_chrono(&param));
	
	return 0;		
}


/*---------------------------------------------------------------------------------------------*/
/*!	\fn	int md1D_standard (struct Param_struct *par, struct Efficacites_struct *eff,struct Noms_fichiers *nomfichier)
 *
 *	\brief	Calcul des efficacités en TE ou TM
 */
/*---------------------------------------------------------------------------------------------*/
int md1D_standard (struct Param_struct *par, struct Efficacites_struct *eff,struct Noms_fichiers *nomfichier)
{

	/* Calcul de la matrice S de la surface */
	matrice_S(par);

	/* Calcul des limites des modes propagatifs */
	md1D_propagativ_limits(par, eff);
		
	/* Amplitude du champ incident */
	md1D_incident_field(par,eff);

	
	/* Amplitude des champs diffractés */
	md1D_amplitudes(par->Ai, par->A0, par->Ah, par->S12, par->S22, par);
/*
fd=fopen("chp_in.bin","w");
fwrite(par->A0, sizeof(complex), 2*N+1, fd);
fclose(fd);
*/
	/* Calcul des efficacités */
	md1D_efficacites(par->Ai, par->A0, par->Ah, par, eff);

	/* Ecriture des résultats dans fichier_results */
	md1D_genere_nom_fichier_results(nomfichier->fichier_results, par);
	md1D_ecrire_results(nomfichier->fichier_results, par, eff);

	return 0;
}


/*---------------------------------------------------------------------------------------------*/
/*!	\fn		int md1D_ellipso (struct Param_struct *par, struct Efficacites_struct *eff,struct Noms_fichiers *nomfichier)
 *
 *	\brief	Calcul des amplitudes, efficacités, déphasages en TE et TM, et du dephasage polarimetrique
 *
 *	\todo	A améliorer, faire plus clean !
 */
/*---------------------------------------------------------------------------------------------*/
int md1D_ellipso (struct Param_struct *par, struct Efficacites_struct *eff,struct Noms_fichiers *nomfichier)
{
	int NeffR, NeffT, n;
	FILE *fp;
	int N = par->N;
	double *effR_s, *effR_p, *effT_s, *effT_p, *deltaR_s, *deltaR_p, *deltaR, *deltaT_s, *deltaT_p, *deltaT;
	complex *A0_s, *A0_p, *Ah_s, *Ah_p; 
	double spec_tan_Psi, spec_cos_delta;

	/* Calcul des limites des modes propagatifs */
	md1D_propagativ_limits(par, eff);
		
	/* Amplitude du champ incident */
	md1D_incident_field(par,eff);

	/* Calculs cas TE */
	par->pola = TE;

	if(par->READ_MAT_S){
		md1D_read_mat_S(par);
	}else{
		matrice_S(par);
	}
/*SaveMatrix2file (par->S22, 2*N+1,2*N+1, "Re","stdout");printf("\n");
*/	/* Calcul des amplitudes */
	md1D_amplitudes(par->Ai, par->A0, par->Ah, par->S12, par->S22, par);
	/* Calcul des efficacités */
	md1D_efficacites(par->Ai, par->A0, par->Ah, par, eff);
	/* Récupération des grandeurs */
	int Nmin_super = eff->Nmin_super, Nmax_super = eff->Nmax_super;
	int Nmin_sub = eff->Nmin_sub, Nmax_sub = eff->Nmax_sub;
	NeffR = Nmax_super-Nmin_super+1;
	NeffT = Nmax_sub-Nmin_sub+1;

	A0_s    = (complex *) malloc(sizeof(complex)*(2*N+1));
	A0_p    = (complex *) malloc(sizeof(complex)*(2*N+1));
	Ah_s    = (complex *) malloc(sizeof(complex)*(2*N+1));
	Ah_p    = (complex *) malloc(sizeof(complex)*(2*N+1));
	effR_s  = (double *) malloc(sizeof(double)*NeffR);
	effR_p  = (double *) malloc(sizeof(double)*NeffR);
	deltaR_s = (double *) malloc(sizeof(double)*(2*N+1));
	deltaR_p = (double *) malloc(sizeof(double)*(2*N+1));
	deltaR   = (double *) malloc(sizeof(double)*(2*N+1));
	effT_s  = (double *) malloc(sizeof(double)*NeffT);
	effT_p  = (double *) malloc(sizeof(double)*NeffT);
	deltaT_s = (double *) malloc(sizeof(double)*(2*N+1));
	deltaT_p = (double *) malloc(sizeof(double)*(2*N+1));
	deltaT   = (double *) malloc(sizeof(double)*(2*N+1));

	CopyDbleTab(effR_s, eff->eff_R, NeffR);   /* Efficacité réfléchie */
	CopyDbleTab(effT_s, eff->eff_T, NeffT);   /* Efficacité transmise */
	for (n=-N; n<=N; n++) {
		A0_s[n+N]= par->A0[n+N];      /* Champ réfléchi */
		Ah_s[n+N]= par->Ah[n+N];      /* Champ transmis */
		deltaR_s[n+N] = carg(par->A0[n+N])*180.0/PI;	/* DeltaR_s */
		deltaT_s[n+N] = carg(par->Ah[n+N])*180.0/PI;	/* DeltaR_s */
	}

	/* Calculs cas TM */
/*	par->clock0 = clock();*/ /* réinitialisation des chronomètres */
/*	time(&(par->time0));
*/	par->pola = TM;
/*	par->matrice_T = matrice_T_TM;
*/	if(par->READ_MAT_S){
		sprintf(par->mat_S_file,"S_TM_%s.txt",par->mat_S_name);
		md1D_read_mat_S(par);
	}else{
		matrice_S(par);
	}
/*	S12_TM = par->S12;
*/	/* Calcul des amplitudes */
	md1D_amplitudes(par->Ai, par->A0, par->Ah, par->S12, par->S22, par);
	/* Calcul des efficacités */
	md1D_efficacites(par->Ai, par->A0, par->Ah, par, eff);
	/* Récupération des grandeurs */
	CopyDbleTab(effR_p, eff->eff_R, NeffR);   /* Efficacité réfléchie */
	CopyDbleTab(effT_p, eff->eff_T, NeffT);   /* Efficacité transmise */
	for (n=-N; n<=N; n++) {
		A0_p[n+N] = par->A0[n+N];      /* Champ réfléchi */
		Ah_p[n+N] = par->Ah[n+N];      /* Champ transmis */
		deltaR_p[n+N] = carg(par->A0[n+N])*180.0/PI;	/* DeltaR_p */
		deltaT_p[n+N] = carg(par->Ah[n+N])*180.0/PI;	/* DeltaR_p */
		deltaR[n+N] = carg(-A0_s[n+N]*conj(A0_p[n+N]))*180.0/PI; 
		deltaT[n+N] = carg(-Ah_s[n+N]*conj(Ah_p[n+N]))*180.0/PI; 
	}

/*	CopyDbleTab(effR_p, eff->eff_R, NeffR);
	for (n=Nmin_super; n<=Nmax_super; n++) {
		A0_p[n-Nmin_super] = par->A0[n+N];
		delta_p[n-Nmin_super] = carg(-par->A0[n+N])*180.0/PI;
		delta[n-Nmin_super] = carg(-A0_s[n-Nmin_super]*conj(A0_p[n-Nmin_super]))*180.0/PI; 
	}
*/
	spec_tan_Psi =	  cabs(A0_p[N]/A0_s[N]);
	spec_cos_delta = -creal(A0_p[N]/A0_s[N])/cabs(A0_p[N]/A0_s[N]);

	/* Ecriture des résultats dans fichier_results */
	md1D_genere_nom_fichier_results(nomfichier->fichier_results, par);
	md1D_ecrire_results(nomfichier->fichier_results, par, eff);

	/* Ajout des résultats ellipsométriques */
	if (!(fp = fopen(nomfichier->fichier_results,"a"))){
		fprintf(stderr, "%s ligne %d : Erreur, impossible d'ouvrir %s\n",__FILE__, __LINE__,nomfichier->fichier_results);
		exit(EXIT_FAILURE);
	}
	int LMAX = 1000; /* NORMALEMENT UNE MACRO */
	
	fprintf(fp,"\neffR_s  = "); ecrire_dble_tab(fp, effR_s,  NeffR, " ", LMAX,"\n");
	fprintf(fp,"\neffR_p  = "); ecrire_dble_tab(fp, effR_p,  NeffR, " ", LMAX,"\n");
	fprintf(fp,"\ndeltaR_s = "); ecrire_dble_tab(fp, deltaR_s, (2*N+1), " ", LMAX,"\n");
	fprintf(fp,"\ndeltaR_p = "); ecrire_dble_tab(fp, deltaR_p, (2*N+1), " ", LMAX,"\n");
	fprintf(fp,"\ndeltaR   = "); ecrire_dble_tab(fp, deltaR,   (2*N+1), " ", LMAX,"\n");
	fprintf(fp,"\neffT_s  = "); ecrire_dble_tab(fp, effT_s,  NeffR, " ", LMAX,"\n");
	fprintf(fp,"\neffT_p  = "); ecrire_dble_tab(fp, effT_p,  NeffR, " ", LMAX,"\n");
	fprintf(fp,"\ndeltaT_s = "); ecrire_dble_tab(fp, deltaT_s, (2*N+1), " ", LMAX,"\n");
	fprintf(fp,"\ndeltaT_p = "); ecrire_dble_tab(fp, deltaT_p, (2*N+1), " ", LMAX,"\n");
	fprintf(fp,"\ndeltaT   = "); ecrire_dble_tab(fp, deltaT,   (2*N+1), " ", LMAX,"\n");
	fprintf(fp,"\nSpecular:");
	fprintf(fp,"\nspec_effR_s = %f",effR_s[-Nmin_super]);
	fprintf(fp,"\nspec_effR_p = %f",effR_p[-Nmin_super]);
	fprintf(fp,"\nspec_deltaR_s = %f",deltaR_s[N]);
	fprintf(fp,"\nspec_deltaR_p = %f",deltaR_p[N]);
	fprintf(fp,"\nspec_deltaR = %f",deltaR[N]);
	fprintf(fp,"\nspec_effT_s = %f",effT_s[-Nmin_sub]);
	fprintf(fp,"\nspec_effT_p = %f",effT_p[-Nmin_sub]);
	fprintf(fp,"\nspec_deltaT_s = %f",deltaT_s[N]);
	fprintf(fp,"\nspec_deltaT_p = %f",deltaT_p[N]);
	fprintf(fp,"\nspec_deltaT = %f",deltaT[N]);

	fprintf(fp,"\nspecR_tan_Psi = %f",spec_tan_Psi);
	fprintf(fp,"\nspecR_cos_delta = %f",spec_cos_delta);
		
	fclose(fp);

	free(effR_s);
	free(effR_p);
	free(deltaR_s);
	free(deltaR_p);
	free(deltaR);
	free(effT_s);
	free(effT_p);
	free(deltaT_s);
	free(deltaT_p);
	free(deltaT);
	free(A0_s);
	free(A0_p);
	free(Ah_s);
	free(Ah_p);

	return 0;
}

/*---------------------------------------------------------------------------------------------*/
/*!	\fn		int md1D_aleat_T_ellipso (struct Param_struct *par, struct Efficacites_struct *eff,struct Noms_fichiers *nomfichier)
 *
 *	\brief	Calcul de diffraction par combinaison partiellement aléatoire de matrices T, 
 * 		les matrices T sont calculées, pour une structure donnée, puis recombinées en utilisant 
 * 		l'algorithme matrices S. Utile pour accéder à des épaisseurs importantes dans le cas
 *		de la diffusion de volume.
 *
 *		Les matrices T sont recombinées par groupes, pour ne pas segmenter trop le volume
 */
/*---------------------------------------------------------------------------------------------*/
int md1D_aleat_T_ellipso (struct Param_struct *par, struct Efficacites_struct *eff,struct Noms_fichiers *nomfichier)
{
	int NeffR, i,j,k, n, n_seg_aleat;
	FILE *fp;
	int N_debut, N = par->N;
	double aleat, h_init, *effR_s, *effR_p, *delta_s, *delta_p, *delta;
	complex  *A0_s, *A0_p;
	
	/* paramètres spécifiques */	
	
	/* vérification que n_sub == n_super (sinon algo matrices T à revoir) */
	if (par->n_sub != par->n_super){
		fprintf(stderr, "ERREUR, %s ligne %d : il faut que n_substrat = n_superstrat dans le cas ALEAT_T_ELLIPSO\n",__FILE__,__LINE__);
		exit(EXIT_FAILURE);
	}
	if (par->h < par->L_segment + par->ecart_type_segment){
		fprintf(stderr, "ERREUR, %s ligne %d : h < L_segment + ecart_type_segment !\n",__FILE__,__LINE__);
		exit(EXIT_FAILURE);
	}
	
	/* Calcul de Nstep_S_total : nb de découpages en matrices T pour obtenir l'épaisseur voulue */
	par->NS_total = ROUND(par->NS*(par->h_total_aleat_T/par->h));
printf("NS       : %d \n",	par->NS);
printf("NS_total : %d \n",	par->NS_total);

	/* Longueur moyenne et ecart type d'un segment en nombre de matrices T */
	int D0_segment = ROUND (par->L_segment*par->NS/par->h);
	int segment_ecart = ROUND (par->ecart_type_segment*par->NS/par->h);
printf("D0_segment : %d \n", D0_segment);
printf("segment_ecart : %d \n", segment_ecart);

	/* Allocations de mémoire */
	int N_segments_estime = par->NS_total/MAX((D0_segment-segment_ecart),1);
printf("N_segments_estime : %d \n", N_segments_estime);
	int *D_segment = malloc(sizeof(int)*N_segments_estime);
	par->sequence_T = malloc(sizeof(int)*par->NS_total);
	
		/*----- Détermination de la séquence des matrices T -----*/
	srand(time(NULL));
	/* Longueurs (en nombre de matrices T) des segments successifs */
	int NS_tmp = 0;
	int N_segments = 0;
	do{
		aleat = 2*((double) rand()/RAND_MAX)-1; /* entre -1 et 1 */
		D_segment[N_segments] = MAX(1, D0_segment + (int)(segment_ecart*aleat));		
		NS_tmp += D_segment[N_segments];
		N_segments++;
	}while(NS_tmp < par->NS_total - D0_segment/2);
	/* Ajustement de la taille à Nstep_S_total */
	while (NS_tmp < par->NS_total){
		/* Ajouter 1 à un segment au hasard */
		n_seg_aleat = ROUND((N_segments+1)*rand()/RAND_MAX);
		D_segment[n_seg_aleat]++;
		NS_tmp++;
	}
	while(NS_tmp > par->NS_total){
		/* Enlever 1 à un segment au hasard */
		n_seg_aleat = ROUND((N_segments+1)*rand()/RAND_MAX);
		if (D_segment[n_seg_aleat] > 1){
			D_segment[n_seg_aleat]--;
			NS_tmp--;
		}
	}
printf("\nD_segments : ");
for(i=0;i<=N_segments-1;i++){
printf("%d ",D_segment[i]);
}
	
printf("\nsequence_T[k] : ");
	/* Assignation de n° de matrices T aux segments */
	k=0;
	for(i=0;i<=N_segments-1;i++){
		aleat = (double) rand()/RAND_MAX; /* entre 0 et 1 */
		N_debut = ROUND(aleat*(par->NS - D_segment[i]));
printf("N_debut = %d \n",	N_debut);
		for(j=0;j<=D_segment[i]-1;j++){
			par->sequence_T[k++] = N_debut + j;
printf("%d ",par->sequence_T[k-1]);
		}
	}

	/* Calcul des limites des modes propagatifs */
	md1D_propagativ_limits(par, eff);
		
	/* Amplitude du champ incident */
	md1D_incident_field(par,eff);
	
	/* Calculs cas TE */
	par->pola = TE;
	matrice_S_aleat_T(par);
/*SaveMatrix2file (par->S22, 2*N+1,2*N+1, "Re","stdout");printf("\n");
*/	/* Calcul des amplitudes */
	md1D_amplitudes(par->Ai, par->A0, par->Ah, par->S12, par->S22, par);
	/* Nouvelle valeur de h */
	h_init = par->h;
	par->h = par->h_total_aleat_T;
	/* Calcul des efficacités */
	md1D_efficacites(par->Ai, par->A0, par->Ah, par, eff);
	/* Valeaur initiale pour h */
	par->h = h_init;
	/* Récupération des grandeurs */
	int Nmin_super = eff->Nmin_super, Nmax_super = eff->Nmax_super;
	NeffR = Nmax_super-Nmin_super+1;

	A0_s    = (complex *) malloc(sizeof(complex)*NeffR);
	A0_p    = (complex *) malloc(sizeof(complex)*NeffR);
	effR_s  = (double *) malloc(sizeof(double)*NeffR);
	effR_p  = (double *) malloc(sizeof(double)*NeffR);
	delta_s = (double *) malloc(sizeof(double)*NeffR);
	delta_p = (double *) malloc(sizeof(double)*NeffR);
	delta   = (double *) malloc(sizeof(double)*NeffR);

	CopyDbleTab(effR_s, eff->eff_R, NeffR);       /* Efficacité réfléchi */
	for (n=Nmin_super; n<=Nmax_super; n++) {
		A0_s[n-Nmin_super] = par->A0[n+N];		/* Champ */
		delta_s[n-Nmin_super] = carg(par->A0[n+N])*180.0/PI;	/* Delta_s */
	}

	/* Calculs cas TM */
/*	par->clock0 = clock();*/ /* réinitialisation des chronomètres */
/*	time(&(par->time0));
*/	par->pola = TM;
/*	par->matrice_T = matrice_T_TM;
*/	matrice_S_aleat_T(par);
/*	S12_TM = par->S12;
*/	/* Calcul des amplitudes */
	md1D_amplitudes(par->Ai, par->A0, par->Ah, par->S12, par->S22, par);
	/* Nouvelle valeur de h */
	par->h = par->h_total_aleat_T;
	/* Calcul des efficacités */
	md1D_efficacites(par->Ai, par->A0, par->Ah, par, eff);
	/* Valeaur initiale pour h */
	par->h = h_init;
	/* Récupération des grandeurs */
	CopyDbleTab(effR_p, eff->eff_R, NeffR);       /* Efficacité */
	for (n=Nmin_super; n<=Nmax_super; n++) {
		A0_p[n-Nmin_super] = par->A0[n+N];		/* Champ */
		delta_p[n-Nmin_super] = carg(par->A0[n+N])*180.0/PI;	/* Delta_p */
		delta[n-Nmin_super] = carg(A0_s[n-Nmin_super]*conj(A0_p[n-Nmin_super]))*180.0/PI; 
	}

	/* Ecriture des résultats dans fichier_results */
	md1D_genere_nom_fichier_results(nomfichier->fichier_results, par);
	md1D_ecrire_results(nomfichier->fichier_results, par, eff);

	/* Ajout des résultats ellipsométriques */
	if (!(fp = fopen(nomfichier->fichier_results,"a"))){
		fprintf(stderr, "%s ligne %d : Erreur, impossible d'ouvrir %s\n",__FILE__, __LINE__,nomfichier->fichier_results);
		exit(EXIT_FAILURE);
	}
	int LMAX = 1000; /* NORMALEMENT UNE MACRO */
	
	fprintf(fp,"\neffR_s  = "); ecrire_dble_tab(fp, effR_s,  NeffR, " ", LMAX,"\n");
	fprintf(fp,"\neffR_p  = "); ecrire_dble_tab(fp, effR_p,  NeffR, " ", LMAX,"\n");
	fprintf(fp,"\ndelta_s = "); ecrire_dble_tab(fp, delta_s, NeffR, " ", LMAX,"\n");
	fprintf(fp,"\ndelta_p = "); ecrire_dble_tab(fp, delta_p, NeffR, " ", LMAX,"\n");
	fprintf(fp,"\ndelta   = "); ecrire_dble_tab(fp, delta,   NeffR, " ", LMAX,"\n");
		
	fclose(fp);

	free(effR_s);
	free(effR_p);
	free(delta_s);
	free(delta_p);
	free(delta);

	free(par->sequence_T);
	free(D_segment);
	
	
	return 0;
}

/*---------------------------------------------------------------------------------------------*/
/*!	\fn		int md1D_replic_T_ellipso (struct Param_struct *par, struct Efficacites_struct *eff,struct Noms_fichiers *nomfichier)
 *
 *	\brief	Calcul de diffraction par replication de sequence de matrice T. Application aux calculs de transmission par la cornee
 *
 */
/*---------------------------------------------------------------------------------------------*/
int md1D_replic_T_ellipso (struct Param_struct *par, struct Efficacites_struct *eff,struct Noms_fichiers *nomfichier)
{
	int NeffR,j,k, n, half_NS;
	FILE *fp;
	int N = par->N, NS=par->NS;
	double h_init, *effR_s, *effR_p, *delta_s, *delta_p, *delta, effR0_s, effR0_p;
	complex  *A0_s, *A0_p;
	
	/* paramètres spécifiques */	
	
	/* vérification que NS est pair */
	if (fabs((double)(par->NS)/2.0-ROUND((double)(par->NS)/2.0))>0.0001){
		fprintf(stderr, "ERROR, %s line %d : NS must be pair for replic_T_ellispo\n",__FILE__,__LINE__);
		exit(EXIT_FAILURE);
	}
	
	/* Calcul de Nstep_S_total : nb de découpages en matrices T pour obtenir l'épaisseur voulue */
	par->NS_total = ROUND(par->NS*(par->h_total_aleat_T/par->h));

printf("NS       : %d \n",	par->NS);
printf("NS_total : %d \n",	par->NS_total);


	/* Allocations de mémoire */
	par->sequence_T = malloc(sizeof(int)*par->NS_total);
	
		/*----- Détermination de la séquence des matrices T -----*/
	half_NS = ROUND(NS/2);
	k=0;
	for(j=NS;j>=half_NS+1;j--){
		par->sequence_T[k++] = j;
	}
	while(k <= par->NS_total-1){
		for(j=half_NS;j>=1 && k <= par->NS_total-1;j--){
			par->sequence_T[k++] = j;
		}
	}
	
	/* Calcul des limites des modes propagatifs */
	md1D_propagativ_limits(par, eff);
		
	/* Amplitude du champ incident */
	md1D_incident_field(par,eff);
	
	/* Calculs cas TE */
	par->pola = TE;
	matrice_S_aleat_T(par);
/*SaveMatrix2file (par->S22, 2*N+1,2*N+1, "Re","stdout");printf("\n");
*/	/* Calcul des amplitudes */
	md1D_amplitudes(par->Ai, par->A0, par->Ah, par->S12, par->S22, par);
	/* Nouvelle valeur de h */
	h_init = par->h;
	par->h = par->h_total_aleat_T;
	/* Calcul des efficacités */
	md1D_efficacites(par->Ai, par->A0, par->Ah, par, eff);
	/* Valeaur initiale pour h */
	par->h = h_init;
	/* Récupération des grandeurs */
	int Nmin_super = eff->Nmin_super, Nmax_super = eff->Nmax_super;
	NeffR = Nmax_super-Nmin_super+1;

	effR0_s = eff->eff_R[-eff->Nmin_super];       /* Efficacité faisceau réfléchi */
	A0_s    = (complex *) malloc(sizeof(complex)*NeffR);
	A0_p    = (complex *) malloc(sizeof(complex)*NeffR);
	effR_s  = (double *) malloc(sizeof(double)*NeffR);
	effR_p  = (double *) malloc(sizeof(double)*NeffR);
	delta_s = (double *) malloc(sizeof(double)*NeffR);
	delta_p = (double *) malloc(sizeof(double)*NeffR);
	delta   = (double *) malloc(sizeof(double)*NeffR);

	CopyDbleTab(effR_s, eff->eff_R, NeffR);       /* Efficacité réfléchi */
	for (n=Nmin_super; n<=Nmax_super; n++) {
		A0_s[n-Nmin_super] = par->A0[n+N];		/* Champ */
		delta_s[n-Nmin_super] = carg(par->A0[n+N])*180.0/PI;	/* Delta_s */
	}

	/* Calculs cas TM */
/*	par->clock0 = clock();*/ /* réinitialisation des chronomètres */
/*	time(&(par->time0));
*/	par->pola = TM;
/*	par->matrice_T = matrice_T_TM;
*/	matrice_S_aleat_T(par);
/*	S12_TM = par->S12;
*/	/* Calcul des amplitudes */
	md1D_amplitudes(par->Ai, par->A0, par->Ah, par->S12, par->S22, par);
	/* Nouvelle valeur de h */
	par->h = par->h_total_aleat_T;
	/* Calcul des efficacités */
	md1D_efficacites(par->Ai, par->A0, par->Ah, par, eff);
	/* Valeaur initiale pour h */
	par->h = h_init;
	/* Récupération des grandeurs */
	CopyDbleTab(effR_p, eff->eff_R, NeffR);       /* Efficacité */
	for (n=Nmin_super; n<=Nmax_super; n++) {
		A0_p[n-Nmin_super] = par->A0[n+N];		/* Champ */
		delta_p[n-Nmin_super] = carg(par->A0[n+N])*180.0/PI;	/* Delta_p */
		delta[n-Nmin_super] = carg(A0_s[n-Nmin_super]*conj(A0_p[n-Nmin_super]))*180.0/PI; 
	}

	/* Ecriture des résultats dans fichier_results */
	md1D_genere_nom_fichier_results(nomfichier->fichier_results, par);
	md1D_ecrire_results(nomfichier->fichier_results, par, eff);

	/* Ajout des résultats ellipsométriques */
	if (!(fp = fopen(nomfichier->fichier_results,"a"))){
		fprintf(stderr, "%s ligne %d : Erreur, impossible d'ouvrir %s\n",__FILE__, __LINE__,nomfichier->fichier_results);
		exit(EXIT_FAILURE);
	}
	int LMAX = 1000; /* NORMALEMENT UNE MACRO */

/*	delta_s[par->ni] = carg(S12_TE[par->N][par->N])*180.0/PI;*/	
	effR0_p = eff->eff_R[-eff->Nmin_super];       /* Efficacité faisceau réfléchi */

	
	fprintf(fp,"\neffR0_s  = %1.12e\n",effR0_s);
	fprintf(fp,"\neffR0_p  = %1.12e\n",effR0_p);
	fprintf(fp,"\neffR_s  = "); ecrire_dble_tab(fp, effR_s,  NeffR, " ", LMAX,"\n");
	fprintf(fp,"\neffR_p  = "); ecrire_dble_tab(fp, effR_p,  NeffR, " ", LMAX,"\n");
	fprintf(fp,"\ndelta_s = "); ecrire_dble_tab(fp, delta_s, NeffR, " ", LMAX,"\n");
	fprintf(fp,"\ndelta_p = "); ecrire_dble_tab(fp, delta_p, NeffR, " ", LMAX,"\n");
	fprintf(fp,"\ndelta   = "); ecrire_dble_tab(fp, delta,   NeffR, " ", LMAX,"\n");
		
	fclose(fp);

	free(effR_s);
	free(effR_p);
	free(delta_s);
	free(delta_p);
	free(delta);

	free(par->sequence_T);
	
	
	return 0;
}

/*---------------------------------------------------------------------------------------------*/
/*!	\fn		int md1D_variation_incidence_A (struct Param_struct *par, struct Efficacites_struct *eff, struct Noms_fichiers *nomfichier)
 *
 *	\brief	Fait varier l'incidence en utilisant une seule matrice S et plusieurs vecteurs chp incident \n
 *			Ai. Adapté aux grands N (ex : calculs de diffusion)
 */
/*---------------------------------------------------------------------------------------------*/
int md1D_variation_incidence_A (struct Param_struct *par, struct Efficacites_struct *eff, struct Noms_fichiers *nomfichier) 
{
	int n,j;
	FILE *fp;
	int N = par->N;
	
	/* Calcul de la matrice S de la surface */
	matrice_S(par);


	/* Appel de md1D_efficacités pour calcul de Nmax_super et Nmin_super */
	md1D_efficacites(par->Ai, par->A0, par->Ah, par, eff);

	for (n=eff->Nmin_super; n<=eff->Nmax_super; n++) {
		for (j=-N; j<=N; j++) {
			par->Ai[j+N] = 0;
		}
		par->Ai[n+N] = 1;
		
		md1D_amplitudes(par->Ai, par->A0, par->Ah, par->S12, par->S22, par);
		md1D_efficacites(par->Ai, par->A0, par->Ah, par, eff);
		
		/* Récupération de la grandeur */
		par->var_i[n-eff->Nmin_super] = eff->eff_R[n-eff->Nmin_super]; /* Efficacité faisceau réfléchi */
		
	}

	/* Ecriture des résultats dans fichier_results */
	md1D_genere_nom_fichier_results(nomfichier->fichier_results, par);
	md1D_ecrire_results(nomfichier->fichier_results, par, eff);

	/* Ajout d'une ligne contenant var_i */
	if (!(fp = fopen(nomfichier->fichier_results,"a"))){
		fprintf(stderr, "%s ligne %d : Erreur, impossible d'ouvrir %s\n",__FILE__, __LINE__,nomfichier->fichier_results);
		exit(EXIT_FAILURE);
	}
	int LMAX = 100; /* NORMALEMENT UNE MACRO */
	fprintf(fp,"\nvar_i = "); ecrire_dble_tab(fp, par->var_i, eff->Nmax_super-eff->Nmin_super+1, " ", LMAX,"\n");
	fclose(fp);


	return 0;
}


/*---------------------------------------------------------------------------------------------*/
/*!	\fn		int md1D_variation_incidence_B (struct Param_struct *par, struct Efficacites_struct *eff, struct Noms_fichiers *nomfichier) 
 *
 *	\brief	Fait varier l'incidence en calculant une matrice S ppour chaque theta_i. \n
 *			Adapté aux petits N (ex : réseaux diélectriques de faibles pas)
 */
/*---------------------------------------------------------------------------------------------*/
int md1D_variation_incidence_B (struct Param_struct *par, struct Efficacites_struct *eff, struct Noms_fichiers *nomfichier) 
{
	FILE *fp;
	double i_min = 0;
	double i_max = 89;
	double delta_i = 1;
	par->Ni = ROUND((i_max - i_min + 1)/delta_i); 
	int Ni = par->Ni;
	double *var_i_angle, *var_i_effR, *var_i_effR_p1, *var_i_effR_m1, *var_i_effR_m2, *var_i_modR, *var_i_argR;

	var_i_angle = (double *) malloc(sizeof(double)*Ni);
	var_i_effR  = (double *) malloc(sizeof(double)*Ni);
	var_i_effR_p1  = (double *) malloc(sizeof(double)*Ni);
	var_i_effR_m1  = (double *) malloc(sizeof(double)*Ni);
	var_i_effR_m2  = (double *) malloc(sizeof(double)*Ni);
	var_i_modR  = (double *) malloc(sizeof(double)*Ni);
	var_i_argR  = (double *) malloc(sizeof(double)*Ni);


	/* Boucle sur l'angle d'incidence */
	for (par->ni=0; par->ni<=Ni-1; (par->ni)++){
	
		par->theta_i = (i_min + delta_i*(double)par->ni)*PI/180.0;
		par->sigma0 = par->k0*sin(par->theta_i);
		var_i_angle[par->ni] = (i_min + delta_i*(double)par->ni);
			
/*		fprintf(stdout,"\r i = %3.0f   ",par->theta_i*180.0/PI);fflush(stdout);
*/		par->verbose = 0;
	
		/* Calcul de la matrice S de la surface */
		matrice_S(par);

		/* Calcul des limites des modes propagatifs */
		md1D_propagativ_limits(par, eff);
			
		/* Amplitude du champ incident */
		md1D_incident_field(par,eff);

		md1D_amplitudes(par->Ai, par->A0, par->Ah, par->S12, par->S22, par);

		/* Calcul des efficacités */
		md1D_efficacites(par->Ai, par->A0, par->Ah, par, eff);
			
		/* Récupération de la grandeur */
		var_i_effR[par->ni] = eff->eff_R[-eff->Nmin_super];       /* Efficacité faisceau réfléchi */
		var_i_effR_p1[par->ni] = eff->eff_R[-eff->Nmin_super+1];       /* Efficacité ordre 1 */
		var_i_effR_m1[par->ni] = eff->eff_R[-eff->Nmin_super-1];       /* Efficacité ordre -1 */
		var_i_effR_m2[par->ni] = eff->eff_R[-eff->Nmin_super-2];       /* Efficacité ordre -1 */
		var_i_modR[par->ni] = cabs(par->S12[par->N][par->N]);     /* module du facteur de réflexion du spéculaire */
		var_i_argR[par->ni] = carg(par->S12[par->N][par->N]); /* argument du facteur de réflexion complexe du spéculaire */

	}
	
	/* Ecriture des résultats dans fichier_results */
	md1D_genere_nom_fichier_results(nomfichier->fichier_results, par);
	md1D_ecrire_results(nomfichier->fichier_results, par, eff);

	/* Ajout d'une ligne contenant var_i */
	if (!(fp = fopen(nomfichier->fichier_results,"a"))){
		fprintf(stderr, "%s ligne %d : Erreur, impossible d'ouvrir %s\n",__FILE__, __LINE__,nomfichier->fichier_results);
		exit(EXIT_FAILURE);
	}
	int LMAX = 1000; /* NORMALEMENT UNE MACRO */
	
	fprintf(fp,"\nvar_i_angle = "); ecrire_dble_tab(fp, var_i_angle, par->Ni, " ", LMAX,"\n");
	fprintf(fp,"\nvar_i_effR  = "); ecrire_dble_tab(fp, var_i_effR, par->Ni, " ", LMAX,"\n");
	fprintf(fp,"\nvar_i_effR_p1  = "); ecrire_dble_tab(fp, var_i_effR_p1, par->Ni, " ", LMAX,"\n");
	fprintf(fp,"\nvar_i_effR_m1  = "); ecrire_dble_tab(fp, var_i_effR_m1, par->Ni, " ", LMAX,"\n");
	fprintf(fp,"\nvar_i_effR_m2  = "); ecrire_dble_tab(fp, var_i_effR_m2, par->Ni, " ", LMAX,"\n");
	fprintf(fp,"\nvar_i_modR  = "); ecrire_dble_tab(fp, var_i_modR, par->Ni, " ", LMAX,"\n");
	fprintf(fp,"\nvar_i_argR  = "); ecrire_dble_tab(fp, var_i_argR, par->Ni, " ", LMAX,"\n");
		
	fclose(fp);


	return 0;
}

/*---------------------------------------------------------------------------------------------*/
/*!	\fn		int md1D_var_i_ellipso (struct Param_struct *par, struct Efficacites_struct *eff,struct Noms_fichiers *nomfichier)
 *
 *	\brief	Calcul des amplitudes, efficacités, déphasages en TE et TM, et du dephasage polarimetrique
 *
 *	\todo	A améliorer, faire plus clean !
 */
/*---------------------------------------------------------------------------------------------*/
int md1D_var_i_ellipso (struct Param_struct *par, struct Efficacites_struct *eff, struct Noms_fichiers *nomfichier)
{
	FILE *fp;
	double i_min = 0;
	double i_max = 89;
	double delta_i = 1;
	par->Ni = ROUND((i_max - i_min +1)/delta_i);
	int Ni = par->Ni;
	double *var_i_angle, *var_i_effR_s, *var_i_effR_p, *var_i_delta_s, *var_i_delta_p, *var_i_delta;

	var_i_angle   = (double *) malloc(sizeof(double)*Ni);
	var_i_effR_s  = (double *) malloc(sizeof(double)*Ni);
	var_i_effR_p  = (double *) malloc(sizeof(double)*Ni);
	var_i_delta_s = (double *) malloc(sizeof(double)*Ni);
	var_i_delta_p = (double *) malloc(sizeof(double)*Ni);
	var_i_delta   = (double *) malloc(sizeof(double)*Ni);

	complex **S12_TE, **S12_TM;
	S12_TE = allocate_CplxMatrix(2*par->N+1, 2*par->N+1);

	/* Calcul des limites des modes propagatifs */
	md1D_propagativ_limits(par, eff);
		
	/* Amplitude du champ incident */
	md1D_incident_field(par,eff);

	/* Boucle sur l'angle d'incidence */
	for (par->ni=0; par->ni<=Ni-1; (par->ni)++){
	
		par->theta_i = (i_min + delta_i*(double)par->ni)*PI/180.0;
		par->sigma0 = par->k0*sin(par->theta_i);
		var_i_angle[par->ni] = (i_min + delta_i*(double)par->ni);
			
/*		fprintf(stdout,"\r i = %3.0f   ",par->theta_i*180.0/PI);fflush(stdout);
*/		par->verbose = 0;

		/* Calculs cas TE */
		par->pola = TE;
/*		par->matrice_T = matrice_T_TE;*/
		matrice_S(par);
		M_egal(S12_TE, par->S12, 2*par->N+1, 2*par->N+1);
		/* Calcul des amplitudes */
		md1D_amplitudes(par->Ai, par->A0, par->Ah, par->S12, par->S22, par);
		/* Calcul des efficacités */
		md1D_efficacites(par->Ai, par->A0, par->Ah, par, eff);
		/* Récupération des grandeurs */
		var_i_effR_s[par->ni] = eff->eff_R[-eff->Nmin_super];       /* Efficacité faisceau réfléchi */
		var_i_delta_s[par->ni] = carg(S12_TE[par->N][par->N])*180.0/PI;       

		/* Calculs cas TM */
		par->pola = TM;
/*		par->matrice_T = matrice_T_TM;
*/		matrice_S(par);
		S12_TM = par->S12;
		/* Calcul des amplitudes */
		md1D_amplitudes(par->Ai, par->A0, par->Ah, par->S12, par->S22, par);
		/* Calcul des efficacités */
		md1D_efficacites(par->Ai, par->A0, par->Ah, par, eff);
		/* Récupération des grandeurs */
		var_i_effR_p[par->ni] = eff->eff_R[-eff->Nmin_super];       /* Efficacité faisceau réfléchi */
		var_i_delta_p[par->ni] = carg(S12_TM[par->N][par->N])*180.0/PI;       
		var_i_delta[par->ni] = carg(S12_TE[par->N][par->N]*conj(S12_TM[par->N][par->N]))*180.0/PI;       
	}

	/* Ecriture des résultats dans fichier_results */
	md1D_genere_nom_fichier_results(nomfichier->fichier_results, par);
	md1D_ecrire_results(nomfichier->fichier_results, par, eff);

	/* Ajout d'une ligne contenant var_i */
	if (!(fp = fopen(nomfichier->fichier_results,"a"))){
		fprintf(stderr, "%s ligne %d : Erreur, impossible d'ouvrir %s\n",__FILE__, __LINE__,nomfichier->fichier_results);
		exit(EXIT_FAILURE);
	}
	int LMAX = 1000; /* NORMALEMENT UNE MACRO */
	
	fprintf(fp,"\nvar_i_angle   = "); ecrire_dble_tab(fp, var_i_angle,   Ni, " ", LMAX,"\n");
	fprintf(fp,"\nvar_i_effR_s  = "); ecrire_dble_tab(fp, var_i_effR_s,  Ni, " ", LMAX,"\n");
	fprintf(fp,"\nvar_i_effR_p  = "); ecrire_dble_tab(fp, var_i_effR_p,  Ni, " ", LMAX,"\n");
	fprintf(fp,"\nvar_i_delta_s = "); ecrire_dble_tab(fp, var_i_delta_s, Ni, " ", LMAX,"\n");
	fprintf(fp,"\nvar_i_delta_p = "); ecrire_dble_tab(fp, var_i_delta_p, Ni, " ", LMAX,"\n");
	fprintf(fp,"\nvar_i_delta   = "); ecrire_dble_tab(fp, var_i_delta,   Ni, " ", LMAX,"\n");
		
	fclose(fp);

	free(S12_TE);

	return 0;
}

/*---------------------------------------------------------------------------------------------*/
/*!	\fn		int md1D_var_lambda_ellipso (struct Param_struct *par, struct Efficacites_struct *eff,struct Noms_fichiers *nomfichier)
 *
 *	\brief	Calcul des amplitudes, efficacités, déphasages en TE et TM, et du dephasage polarimetrique
 *
 *	\todo	A améliorer, faire plus clean !
 */
/*---------------------------------------------------------------------------------------------*/
int md1D_var_lambda_ellipso (struct Param_struct *par, struct Efficacites_struct *eff, struct Noms_fichiers *nomfichier)
{
	int i;
	FILE *fp;
	int N_layers = par->N_layers;
	char material[10][SIZE_STR_BUFFER];
	double lambda_min = 210;
	double lambda_max = 780;
	double delta_lambda = 10;
	par->Ni = ROUND((lambda_max - lambda_min +1)/delta_lambda);
	int Ni = par->Ni;
	double *var_lambda, *var_tan_Psi, *var_cos_delta, *var_i_effR_s, *var_i_effR_p, *var_i_delta_s, *var_i_delta_p, *var_i_delta;
	complex indice;

	if (par->type_profil != H_X && par->type_profil != MULTICOUCHES) {
		fprintf(stderr,"type de profil incompatible avec fonction var_lambda");
		return -1;}
		
	sprintf(material[0],"AIR");
	sprintf(material[1],"POLY03");
	sprintf(material[2],"OXIDE_THERM");
	sprintf(material[3],"SI_CRISTAL");
	
	var_lambda    = (double *) malloc(sizeof(double)*Ni);
	var_tan_Psi   = (double *) malloc(sizeof(double)*Ni);
	var_cos_delta = (double *) malloc(sizeof(double)*Ni);
	var_i_effR_s  = (double *) malloc(sizeof(double)*Ni);
	var_i_effR_p  = (double *) malloc(sizeof(double)*Ni);
	var_i_delta_s = (double *) malloc(sizeof(double)*Ni);
	var_i_delta_p = (double *) malloc(sizeof(double)*Ni);
	var_i_delta   = (double *) malloc(sizeof(double)*Ni);

	complex **S12_TE, **S12_TM;
	S12_TE = allocate_CplxMatrix(2*par->N+1, 2*par->N+1);

	/* Calcul des limites des modes propagatifs */
	md1D_propagativ_limits(par, eff);
		
	/* Amplitude du champ incident */
	md1D_incident_field(par,eff);

	/* Boucle sur lambda */
	for (par->ni=0; par->ni<=Ni-1; (par->ni)++){
	
		par->lambda = lambda_min + delta_lambda*(double)par->ni;
		par->n_sub = md1D_indice("SI_CRISTAL",par->lambda,"Lookup");
		var_lambda[par->ni] = par->lambda;
		par->k0 = 2*PI*par->n_super/par->lambda;
		par->kh = 2*PI*par->n_sub/par->lambda;
		par->sigma0 = par->k0*sin(par->theta_i);
		
		/* Détermination des indices pour le lambda considéré */
		for (i=0; i<=N_layers+1; i++){
			indice = md1D_indice(material[i],par->lambda,"Lookup");
			par->k2_layer[i]    = (indice*par->k0)*(indice*par->k0); 
			par->invk2_layer[i] = 1/par->k2_layer[i]; 
		}
		/*Réinitialisation de z2n */
		par->N_z2n = 0;
		for (i=0;i<=par->BLOC_TAILLE_z2n-1;i++){
			par->z2n[0][i] = -939498987;
			par->z2n[1][i] = -998621328;
		}

		
		fprintf(stdout," lambda = %3.0f   ",par->lambda);fflush(stdout);
		par->verbose = 0;

		/* Calculs cas TE */
		par->pola = TE;
/*		par->matrice_T = matrice_T_TE;*/
		matrice_S(par);
		M_egal(S12_TE, par->S12, 2*par->N+1, 2*par->N+1);
		/* Calcul des amplitudes */
		md1D_amplitudes(par->Ai, par->A0, par->Ah, par->S12, par->S22, par);
		/* Calcul des efficacités */
		md1D_efficacites(par->Ai, par->A0, par->Ah, par, eff);
		/* Récupération des grandeurs */
		var_i_effR_s[par->ni] = eff->eff_R[-eff->Nmin_super];       /* Efficacité faisceau réfléchi */
		var_i_delta_s[par->ni] = carg(S12_TE[par->N][par->N])*180.0/PI;       
		
		/* Calculs cas TM */
		par->pola = TM;
/*		par->matrice_T = matrice_T_TM;
*/		matrice_S(par);
		S12_TM = par->S12;
		/* Calcul des amplitudes */
		md1D_amplitudes(par->Ai, par->A0, par->Ah, par->S12, par->S22, par);
		/* Calcul des efficacités */
		md1D_efficacites(par->Ai, par->A0, par->Ah, par, eff);
		/* Récupération des grandeurs */
		var_i_effR_p[par->ni] = eff->eff_R[-eff->Nmin_super];       /* Efficacité faisceau réfléchi */
		var_i_delta_p[par->ni] = carg(S12_TM[par->N][par->N])*180.0/PI;       
		var_i_delta[par->ni] = carg(S12_TE[par->N][par->N]*conj(S12_TM[par->N][par->N]))*180.0/PI;       
		var_tan_Psi[par->ni] = cabs( S12_TM[par->N][par->N] / S12_TE[par->N][par->N] );
		var_cos_delta[par->ni] = creal( S12_TM[par->N][par->N] / S12_TE[par->N][par->N] ) / var_tan_Psi[par->ni];
	}

	/* Ecriture des résultats dans fichier_results */
	md1D_genere_nom_fichier_results(nomfichier->fichier_results, par);
	md1D_ecrire_results(nomfichier->fichier_results, par, eff);

	/* Ajout d'une ligne contenant var_i */
	if (!(fp = fopen(nomfichier->fichier_results,"a"))){
		fprintf(stderr, "%s ligne %d : Erreur, impossible d'ouvrir %s\n",__FILE__, __LINE__,nomfichier->fichier_results);
		exit(EXIT_FAILURE);
	}
	int LMAX = 1000; /* NORMALEMENT UNE MACRO */
	
	fprintf(fp,"\nvar_lambda    = "); ecrire_dble_tab(fp, var_lambda,    Ni, " ", LMAX,"\n");
	fprintf(fp,"\nvar_tan_Psi   = "); ecrire_dble_tab(fp, var_tan_Psi,   Ni, " ", LMAX,"\n");
	fprintf(fp,"\nvar_cos_delta = "); ecrire_dble_tab(fp, var_cos_delta, Ni, " ", LMAX,"\n");
	fprintf(fp,"\nvar_i_effR_s  = "); ecrire_dble_tab(fp, var_i_effR_s,  Ni, " ", LMAX,"\n");
	fprintf(fp,"\nvar_i_effR_p  = "); ecrire_dble_tab(fp, var_i_effR_p,  Ni, " ", LMAX,"\n");
	fprintf(fp,"\nvar_i_delta_s = "); ecrire_dble_tab(fp, var_i_delta_s, Ni, " ", LMAX,"\n");
	fprintf(fp,"\nvar_i_delta_p = "); ecrire_dble_tab(fp, var_i_delta_p, Ni, " ", LMAX,"\n");
	fprintf(fp,"\nvar_i_delta   = "); ecrire_dble_tab(fp, var_i_delta,   Ni, " ", LMAX,"\n");
	fprintf(fp,"\nNb_lambda = %d",par->Ni);
		
	fclose(fp);

	free(var_lambda);
	free(var_tan_Psi);
	free(var_cos_delta);
	free(var_i_effR_s);
	free(var_i_effR_p);
	free(var_i_delta_s);
	free(var_i_delta_p);
	free(var_i_delta);
	free(S12_TE);

	return 0;
}

/*---------------------------------------------------------------------------------------------*/
/*!	\fn	int md1D_standard (struct Param_struct *par, struct Efficacites_struct *eff,struct Noms_fichiers *nomfichier)
 *
 *	\brief	Calcul des efficacités en TE ou TM
 */
/*---------------------------------------------------------------------------------------------*/
int md1D_near_field (struct Param_struct *par, struct Efficacites_struct *eff,struct Noms_fichiers *nomfichier)
{
#if 0
	par->Near_field_Matrix = allocate_CplxMatrix(NS(+1?), 2*par->N+1);
	
	/* 1st iteration : Calculation of the bottom field */
	par->SAVE_NEAR_FIELD = 0;
	matrice_S(par);
	md1D_propagativ_limits(par, eff);
	md1D_incident_field(par,eff);
	md1D_amplitudes(par->Ai, par->A0, par->Ah, par->S12, par->S22, par);
	par->bottom_field = par->Ah;
	
	/* 2nd iteration with field components extraction at each step */
	par->SAVE_NEAR_FIELD = 1;
	matrice_S(par);

	
	/* Calcul des efficacités */
	md1D_efficacites(par->Ai, par->A0, par->Ah, par, eff);

	/* Ecriture des résultats dans fichier_results */
	md1D_genere_nom_fichier_results(nomfichier->fichier_results, par);
	md1D_ecrire_results(nomfichier->fichier_results, par, eff);
#endif
	return 0;
}


/*---------------------------------------------------------------------------------------------*/
/*!	\fn		int md1D_init (struct Param_struct *par, struct Efficacites_struct *eff,struct Noms_fichiers *nomfichier) 
 *
 *	\brief	Initialisation du programme
 */
/*---------------------------------------------------------------------------------------------*/
int md1D_init (struct Param_struct *par, struct Efficacites_struct *eff,struct Noms_fichiers *nomfichier) 
{

	/* Lecture des paramètres par defaut dans fichier_param */
	md1D_lire_param(nomfichier, par);

	/* Allocation de mémoire pour le profil */
	md1D_alloc_init_profil(par);

	/* Lecture du profil h(x) décrivant la surface */
	(*par->md1D_lire_profil)(nomfichier->fichier_profil, par);

	/* Initialisations de certaines variables */
	md1D_variables_init(par, eff);

	/* Allocation de mémoire pour les tableaux */
	md1D_alloc(par, eff);

	/* Affichage des paramètres lus et calculés */
	md1D_affiche_valeurs_param(par, nomfichier);

	return 0;
}



/*---------------------------------------------------------------------------------------------*/
/*!	\fn		int md1D_variables_init(struct Param_struct *par, struct Efficacites_struct *eff)
 *
 *	\brief	Initialisations des variables
 */
/*---------------------------------------------------------------------------------------------*/
int md1D_variables_init(struct Param_struct *par, struct Efficacites_struct *eff)
{
	double K,sig0;
	par->clock0 = clock(); /* Initialisation des chronomètres                  */
	time(&(par->time0));   /* clock0 (courte durées) et time0 (longues durées) */
	par->last_clock = clock();
	time(&(par->last_time)); 

	/* Initialisation des variables en mode "AUTO" */
	/* delta_h = lambda x (n_re + n_im) / 1000 */
	if (par->delta_h == AUTO) {
		par->delta_h = par->lambda/(cabs(par->n_sub)*1000);
	}
	/* NS = (h/lambda) x (n_re + n_im) x 5 */
	if (par->NS == AUTO) {
		par->NS = CEIL( (par->h / par->lambda)*cabs(par->n_sub)*5 );
	}
	/* N : on ajoute 10% de modes evanescents (3 au minimum) */
	if (par->N == AUTO) {
		K = 2.0*PI/par->L;
		if (!strcmp(par->calcul_type,"STD") || !strcmp(par->calcul_type,"ELLIPSO") ){
			sig0 = par->sigma0;
		}else if (!strcmp(par->calcul_type,"VAR_I") || !strcmp(par->calcul_type,"VAR_I_BIS") || \
			!strcmp(par->calcul_type,"VAR_I_ELLIPSO")){
			sig0 = par->k0; /* correspond à sigma0 pour theta_i = 90° => N constant et suffisant de 0 à 90°*/
		}else{
			fprintf(stderr, "%s ligne %d, Calcul AUTO de N : %s, type calcul inconnu\n",__FILE__, __LINE__,par->calcul_type);
			exit(EXIT_FAILURE);
		}
		int N_limit =  FLOOR(( MAX(creal(par->k0),creal(par->kh)) + fabs(sig0))/K);
		/* On ajoute 10% de modes evanescents (3 au minimum) */
		int N_evanesc = ROUND(MAX(3,0.1*N_limit));
		par->N = N_limit + N_evanesc; 
	}

	
	
	
	/* ni et Ni, pour compteurs en theta_i */
	par->ni = 0;
	par->Ni = 1; /* Pour compatibilité (une autre valeur sera affectée par les fonctions var_i_... )*/
	
	/* Calcul de la valeur exacte de delta_h de sorte qu'il y en ait un nb entier à chaque étape Matrice-S*/
	par->Nstep_S = ROUND(ceil((par->h/par->NS)/par->delta_h));
	par->delta_h = (par->h/par->NS)/par->Nstep_S;
	par->Nstep = ROUND(par->h/par->delta_h);

	par->verbose = 2; /* Si VERBOSE > 0 : affiche plus d'infos sur le terminal */
	
	par->STOCKER_TF = 1; /* Optimisation : On stocke ou non les TFs, + rapide, mais demande + de mémoire */
	par->BLOC_TAILLE_z2n = (int) 10*par->Nstep; /* PAS OPTIMISÉ, à REVOIR */
	par->TAILLE_z2n = par->BLOC_TAILLE_z2n;
	par->N_z2n = 0;

	return 0;
}



/*---------------------------------------------------------------------------------------------*/
/*!	\fn	int md1D_alloc_init_profil(struct Param_struct *par)
 *
 *	\brief	Memory allocation for profil
 */
/*---------------------------------------------------------------------------------------------*/
int md1D_alloc_init_profil(struct Param_struct *par)
{
	/* Initilisation de variables */
	par->k0 = 2*PI*par->n_super/par->lambda;
	par->kh = 2*PI*par->n_sub/par->lambda;
	par->delta_sigma = 2*PI/par->L;
	par->sigma0 = par->k0*sin(par->theta_i);

	/* Allocations */	
	if (par->type_profil == N_XYZ) {
		par->n_xyz = allocate_CplxMatrix(par->N_z,par->N_x);
	}else{
		par->profil = allocate_DbleMatrix(par->N_layers+3,par->N_x);
		par->k2_layer    = (complex *) malloc(sizeof(complex)*(par->N_layers+2));
		par->invk2_layer = (complex *) malloc(sizeof(complex)*(par->N_layers+2));
	}

	/* Alignement des pointeurs de fonction */
/*	par->matrice_T = (par->pola == TE ? matrice_T_TE : matrice_T_TM);
*/	switch (par->type_profil) {
		case H_X          : 
			par->md1D_lire_profil = md1D_lire_profil_H_X;
			par->k_2 = k2_H_X;
			par->invk_2 = invk2_H_X;
			break;
		case MULTICOUCHES : 
			par->md1D_lire_profil = md1D_lire_profil_MULTI;
			par->k_2 = k2_MULTI;
			par->invk_2 = invk2_MULTI;
			break;
		case N_XYZ        : 
			par->md1D_lire_profil = md1D_lire_profil_N_XYZ;
			par->k_2 = k2_N_XYZ;
			par->invk_2 = invk2_N_XYZ;
			break;
	}

	
	return 0;
}

/*---------------------------------------------------------------------------------------------*/
/*!	\fn		int md1D_alloc(struct Param_struct *par, struct Efficacites_struct *eff)
 *
 *	\brief	Memory allocation
 */
/*---------------------------------------------------------------------------------------------*/
int md1D_alloc(struct Param_struct *par, struct Efficacites_struct *eff)
{

	par->k2        = (complex *) malloc(sizeof(complex)*par->N_x);
	par->invk2     = (complex *) malloc(sizeof(complex)*par->N_x);
	par->TF_k2     = (complex *) malloc(sizeof(complex)*(4*par->N+1));
	par->TF_invk2  = (complex *) malloc(sizeof(complex)*(4*par->N+1));
	par->S12 = allocate_CplxMatrix(2*par->N+1,2*par->N+1);
	par->S22 = allocate_CplxMatrix(2*par->N+1,2*par->N+1);
	par->Ai = (complex *) malloc(sizeof(complex)*(2*par->N+1));
	par->A0 = (complex *) malloc(sizeof(complex)*(2*par->N+1));
	par->Ah = (complex *) malloc(sizeof(complex)*(2*par->N+1));	

	eff->eff_R =       (double *) malloc(sizeof(double)*(2*par->N+1));
	eff->eff_T =       (double *) malloc(sizeof(double)*(2*par->N+1));
	eff->N_eff_R =     (double *) malloc(sizeof(double)*(2*par->N+1));
	eff->N_eff_T =     (double *) malloc(sizeof(double)*(2*par->N+1));
	eff->theta_eff_R = (double *) malloc(sizeof(double)*(2*par->N+1));
	eff->theta_eff_T = (double *) malloc(sizeof(double)*(2*par->N+1));

	par->var_i =       (double *) malloc(sizeof(double)*1000); /* A REVOIR */
	par->var_i2 =       (double *) malloc(sizeof(double)*1000); /* A REVOIR */

	if (par->STOCKER_TF) {
		par->z2n[0] = (long int *) malloc(sizeof(long int)*(par->BLOC_TAILLE_z2n));
		par->z2n[1] = (long int *) malloc(sizeof(long int)*(par->BLOC_TAILLE_z2n));
		par->tab_TF_k2 = allocate_CplxMatrix(par->BLOC_TAILLE_z2n, 4*par->N+1);
		par->tab_TF_invk2 = allocate_CplxMatrix(par->BLOC_TAILLE_z2n, 4*par->N+1);
	}
	
	return 0;
}



/*---------------------------------------------------------------------------------------------*/
/*!	\fn		int md1D_free(struct Param_struct *par, struct Efficacites_struct *eff)
 *
 *	\brief	Libération de la mémoire
 */
/*---------------------------------------------------------------------------------------------*/
int md1D_free(struct Param_struct *par, struct Efficacites_struct *eff)
{
	if (par->verbose>=3) fprintf(stdout,"Libération de la mémoire : "); fflush(stdout);

	
	if (par->type_profil == N_XYZ) {
		free(par->n_xyz[0]);
		free(par->n_xyz);
	}else{
		free(par->profil[0]);
		free(par->profil);
		free(par->k2_layer);
		free(par->invk2_layer);
	}

	if (par->STOCKER_TF) {
		free(par->z2n[0]);
		free(par->z2n[1]);
		free(par->tab_TF_k2[0]);
		free(par->tab_TF_k2);
		free(par->tab_TF_invk2[0]);
		free(par->tab_TF_invk2);
	}

	free(par->k2);
	free(par->invk2);
	free(par->TF_k2);
	free(par->TF_invk2);
	free(par->S12[0]);
	free(par->S12);
	free(par->S22[0]);
	free(par->S22);
	free(par->A0);
	free(par->Ai);
	free(par->Ah);
	free(eff->eff_R);
	free(eff->eff_T);
	free(eff->N_eff_R);
	free(eff->N_eff_T);
	free(eff->theta_eff_R);
	free(eff->theta_eff_T);
	
	if (par->verbose>=3) fprintf(stdout,"OK\n"); fflush(stdout);

	
	return 0;
}

/*---------------------------------------------------------------------------------------------*/
/*!	\fn	int md1D_read_mat_S(struct Param_struct *par)
 *
 *	\brief	Lecture de matrices S
 */
/*---------------------------------------------------------------------------------------------*/
int md1D_read_mat_S(struct Param_struct *par)
{
	double **re_tmp, **im_tmp;
	int i,j;
	char *nom_fichier = par->mat_S_file;
	char *erreur="NO_ERROR                     ";
	FILE *fp;
	
	/* Ouverture du fichier */
	if (!(fp = fopen(nom_fichier,"r"))){
		fprintf(stderr, "%s ligne %d : ERREUR, impossible d'ouvrir %s\n",__FILE__, __LINE__,nom_fichier);
		return 1;
	}
	/* Lecture des paramètres */
	if (lire_string (fp, "mat_S_name", par->mat_S_name)) erreur="mat_S_name";
	if (lire_double (fp, "h_partial", &(par->h) )) erreur="h_partial";
/*	if (lire_double (fp, "L", &(par->L) )) erreur="L";
	if (lire_double (fp, "Re_k0", &Re_k0 )) erreur="Re_k0";
	if (lire_double (fp, "Re_kh", &Re_kh )) erreur="Re_kh";
	if (lire_double (fp, "Im_k0", &Re_kh )) erreur="Im_k0";
	if (lire_double (fp, "Im_kh", &Re_kh )) erreur="Im_kh";
	if (lire_double (fp, "delta_sigma", &(par->delta_sigma) )) erreur="delta_sigma";
	if (lire_double (fp, "sigma0", &(par->sigma0) )) erreur="sigma0";
	if (lire_int (fp, "N", &(par->N) )) erreur="N";
*/	fclose(fp);
	/* Vérification de l'absence d'erreurs de lecture */
	if (strcmp(erreur,"NO_ERROR                     ")){
		fprintf(stderr, "%s : Erreur, probleme de lecture de \"%s\"\n",__FILE__,erreur);
		exit(EXIT_FAILURE);
	}
	
	/* Lecture des éléments de matrice S */	
	re_tmp = allocate_DbleMatrix(2*par->N+1,2*par->N+1);
	im_tmp = allocate_DbleMatrix(2*par->N+1,2*par->N+1);
	/* S12 */
	lire_tab(nom_fichier, "Re_S12", re_tmp[0], (2*par->N+1)*(2*par->N+1));
	lire_tab(nom_fichier, "Im_S12", im_tmp[0], (2*par->N+1)*(2*par->N+1));
	for (j=0;j<=2*par->N;j++){
		for (i=0;i<=2*par->N;i++){
			par->S12[i][j] = c_omplex(re_tmp[i][j],im_tmp[i][j]);
		}
	}
	/* S22 */
	lire_tab(nom_fichier, "Re_S22", re_tmp[0], (2*par->N+1)*(2*par->N+1));
	lire_tab(nom_fichier, "Im_S22", im_tmp[0], (2*par->N+1)*(2*par->N+1));
	for (j=0;j<=2*par->N;j++){
		for (i=0;i<=2*par->N;i++){
			par->S22[i][j] = c_omplex(re_tmp[i][j],im_tmp[i][j]);
		}
	}

	free(re_tmp[0]);free(re_tmp);
	free(im_tmp[0]);free(im_tmp);
	
	return 0;		
}




