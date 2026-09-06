/*---------------------------------------------------------------------------------------------*/
/*!	\fn		int md2D_comb_mat_S_ellipso (struct Param_struct *par, struct Efficacites_struct *eff,struct Noms_fichiers *nomfichier)
 *
 *	\brief	Calcul des amplitudes, efficacités, déphasages en TE et TM, et du dephasage polarimetrique
 *
 *	\todo	TIRÉ DE md2D_ELLIPSO A améliorer, pas propre, "provisoire" !
 */
/*---------------------------------------------------------------------------------------------*/
int md2D_comb_mat_S_ellipso (struct Param_struct *par, struct Efficacites_struct *eff,struct Noms_fichiers *nomfichier)
{
	int NeffR,i,k,n;
	FILE *fp;
	int N = par->N;
	double *A0_s, *A0_p, *effR_s, *effR_p, *delta_s, *delta_p, *delta;


	complex **S11, **S12, **S21, **S22;
	complex **S11_1, **S12_1, **S21_1, **S22_1;
	complex **S11_2, **S12_2, **S21_2, **S22_2;
	
	par->S11 = allocate_CplxMatrix(2*par->N+1, 2*par->N+1);
	par->S21 = allocate_CplxMatrix(2*par->N+1, 2*par->N+1);

	S11 = par->S11;
	S12 = par->S12;
	S21 = par->S21;
	S22 = par->S22;
	
	S11_1 = allocate_CplxMatrix(2*par->N+1, 2*par->N+1);
	S12_1 = allocate_CplxMatrix(2*par->N+1, 2*par->N+1);
	S21_1 = allocate_CplxMatrix(2*par->N+1, 2*par->N+1);
	S22_1 = allocate_CplxMatrix(2*par->N+1, 2*par->N+1);

	S11_2 = allocate_CplxMatrix(2*par->N+1, 2*par->N+1);
	S12_2 = allocate_CplxMatrix(2*par->N+1, 2*par->N+1);
	S21_2 = allocate_CplxMatrix(2*par->N+1, 2*par->N+1);
	S22_2 = allocate_CplxMatrix(2*par->N+1, 2*par->N+1);
	
	
	/* Champ incident */
	for (k=-N; k<=N; k++) { 
		par->Ai[k+N] = 0;
	}
	par->Ai[N] = 1;

	/* Calculs cas TE */
	par->pola = TE;

	/* Calcul de la matrice S initiale */
	matrice_S_totale(par);
	
	/* Combinaison des matrices S */
	double h_base = par->h; /* Mise en mémoire de h pour futurs calculs TM */
int N_comb_mat_S=1;
	for (i=0;i<=N_comb_mat_S-1;i++){
		/* S_1 = S */ 
		M_egal(S11_1, par->S11, 2*par->N+1, 2*par->N+1);
		M_egal(S12_1, par->S12, 2*par->N+1, 2*par->N+1);
		M_egal(S21_1, par->S21, 2*par->N+1, 2*par->N+1);
		M_egal(S22_1, par->S22, 2*par->N+1, 2*par->N+1);
	
		/* S_2 = S */ 
		M_egal(S11_2, par->S11, 2*par->N+1, 2*par->N+1);
		M_egal(S12_2, par->S12, 2*par->N+1, 2*par->N+1);
		M_egal(S21_2, par->S21, 2*par->N+1, 2*par->N+1);
		M_egal(S22_2, par->S22, 2*par->N+1, 2*par->N+1);
		
		/* S = comb(S1,S2) */
		md2D_comb_mat_S(S11,S12,S21,S22, S11_1,S12_1,S21_1,S22_1, S11_2,S12_2,S21_2,S22_2, 2*N+1);
	
		/* h = h x 2 */
		par->h *= 2;	
	}
		
	/* Calcul des amplitudes */
	md2D_amplitudes(par->Ai, par->A0, par->Ah, par->S12, par->S22, par);
	/* Calcul des efficacités */
	md2D_efficacites(par->Ai, par->A0, par->Ah, par, eff);
	/* Récupération des grandeurs */
	int Nmin_super = eff->Nmin_super, Nmax_super = eff->Nmax_super;
	NeffR = Nmax_super-Nmin_super+1;

	A0_s    = (double *) malloc(sizeof(double)*NeffR);
	A0_p    = (double *) malloc(sizeof(double)*NeffR);
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
	par->h = h_base; /* réinitialisation de h */
	par->pola = TM;

	/* Calcul de la matrice S initiale */
	matrice_S_totale(par);

	/* Combinaison des matrices S */
	for (i=0;i<=N_comb_mat_S-1;i++){
		/* S_1 = S */ 
		M_egal(S11_1, par->S11, 2*par->N+1, 2*par->N+1);
		M_egal(S12_1, par->S12, 2*par->N+1, 2*par->N+1);
		M_egal(S21_1, par->S21, 2*par->N+1, 2*par->N+1);
		M_egal(S22_1, par->S22, 2*par->N+1, 2*par->N+1);
	
		/* S_2 = S */ 
		M_egal(S11_2, par->S11, 2*par->N+1, 2*par->N+1);
		M_egal(S12_2, par->S12, 2*par->N+1, 2*par->N+1);
		M_egal(S21_2, par->S21, 2*par->N+1, 2*par->N+1);
		M_egal(S22_2, par->S22, 2*par->N+1, 2*par->N+1);
		
		/* S = comb(S1,S2) */
		md2D_comb_mat_S(S11,S12,S21,S22, S11_1,S12_1,S21_1,S22_1, S11_2,S12_2,S21_2,S22_2, 2*N+1);
	
		/* h = h x 2 */
		par->h *=2;	
	}

	/* Calcul des amplitudes */
	md2D_amplitudes(par->Ai, par->A0, par->Ah, par->S12, par->S22, par);
	/* Calcul des efficacités */
	md2D_efficacites(par->Ai, par->A0, par->Ah, par, eff);
	/* Récupération des grandeurs */
	CopyDbleTab(effR_p, eff->eff_R, NeffR);       /* Efficacité */
	for (n=Nmin_super; n<=Nmax_super; n++) {
		A0_p[n-Nmin_super] = par->A0[n+N];		/* Champ */
		delta_p[n-Nmin_super] = carg(par->A0[n+N])*180.0/PI;	/* Delta_p */
		delta[n-Nmin_super] = carg(A0_s[n-Nmin_super]*conj(A0_p[n-Nmin_super]))*180.0/PI; 
	}

	/* Ecriture des résultats dans fichier_results */
	md2D_genere_nom_fichier_results(nomfichier->fichier_results, par);
	md2D_ecrire_results(nomfichier->fichier_results, par, eff);

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

	free(S11[0]);
	free(S11);
	free(S21[0]);
	free(S21);
	
	free(S11_1[0]);
	free(S11_1);
	free(S12_1[0]);
	free(S12_1);
	free(S21_1[0]);
	free(S21_1);
	free(S22_1[0]);
	free(S22_1);
	
	free(S11_2[0]);
	free(S11_2);
	free(S12_2[0]);
	free(S12_2);
	free(S21_2[0]);
	free(S21_2);
	free(S22_2[0]);
	free(S22_2);
	
	
	free(effR_s);
	free(effR_p);
	free(delta_s);
	free(delta_p);
	free(delta);

	return 0;
}


/*-------------------------------------------------------------------------------------*/
/*!	\fn		int matrice_S_totale(struct Param_struct *par)
 *
 *	\brief	Calcul de tous les éléments de la matrice S
 *
 *	Utile pour le mode COMB_MAT_S, sinon seuls 2 sur 4 des elts sont calculés
 */
/*-------------------------------------------------------------------------------------*/
int matrice_S_totale(struct Param_struct *par){

	if (par->verbose) fprintf(stdout,"Calcul de la matrice S\n");

	int N = par->N;
	int NS = par->NS;
	complex **S12, **S22, **S11, **S21; 
	S12 = par->S12;
	S22 = par->S22;
	S11 = par->S11;
	S21 = par->S21;

	int i, j, nS;
	complex **T11, **T12, **T21, **T22, **Z, **tmp, **tmp2, **tmp3;
	complex *F_plus, *F_moins, *F_plus2, *F_moins2;
	
	/* Allocations */
	T11 = allocate_CplxMatrix(2*N+1, 2*N+1);
	T12 = allocate_CplxMatrix(2*N+1, 2*N+1);
	T21 = allocate_CplxMatrix(2*N+1, 2*N+1);
	T22 = allocate_CplxMatrix(2*N+1, 2*N+1);
	Z   = allocate_CplxMatrix(2*N+1, 2*N+1);
	tmp = allocate_CplxMatrix(2*N+1, 2*N+1);
	tmp2= allocate_CplxMatrix(2*N+1, 2*N+1);
	tmp3= allocate_CplxMatrix(2*N+1, 2*N+1);
	F_plus   = (complex *) malloc(sizeof(complex)*(4*N+2));
	F_moins  = (complex *) malloc(sizeof(complex)*(4*N+2));
	F_plus2  = (complex *) malloc(sizeof(complex)*(4*N+2));
	F_moins2 = (complex *) malloc(sizeof(complex)*(4*N+2));
	
	/* Initialisations */
	for (i=-N; i<=N; i++) {
		for (j=-N; j<=N; j++) {
			S11[i+N][j+N] = (i==j);
			S12[i+N][j+N] = 0;
			S21[i+N][j+N] = 0;
			S22[i+N][j+N] = (i==j);
		}
	}
	
	/* Itérations */
	for (nS=NS; nS>=1; nS--) {

		/* Calcul de la matrice T */
		if (par->pola == TE){
			matrice_T_TE(T11, T12, T21, T22, F_plus, F_moins, F_plus2, F_moins2, nS, par);
		}else{
			matrice_T_TM(T11, T12, T21, T22, F_plus, F_moins, F_plus2, F_moins2, nS, par);
		}

		/* Z = inv(T11 + T12*S12) */
		invM(Z, add_M(tmp2, 
			T11, M_x_M(tmp,
				T12,S12,	2*N+1,2*N+1),2*N+1,2*N+1),2*N+1);

		/* S12 = (T21 +T22*S12)*Z */
		M_x_M(S12,
			add_M(tmp2, T21, M_x_M(tmp,
					T22,S12,2*N+1,2*N+1),2*N+1,2*N+1),
			Z,2*N+1, 2*N+1);

		/* S22 = S22*Z */
		M_egal(tmp,S22,2*N+1,2*N+1);
		M_x_M(S22,tmp,Z,2*N+1,2*N+1);

		/* S21 = S21 - S22*T12*S11 */
		M_egal(tmp,S21,2*N+1,2*N+1);
		sub_M(S21,
			tmp,
			M_x_M(tmp2,
				S22,
				M_x_M(tmp3,
					T12,
					S11,	2*N+1,2*N+1),2*N+1,2*N+1),2*N+1,2*N+1);

		/* S11 = (T22 - S12*T12)*S11 */
		M_egal(tmp,S11,2*N+1,2*N+1);
		M_x_M(S11,
			sub_M(tmp2,
				T22,
				M_x_M(tmp3,
					S12,
					T12,	2*N+1,2*N+1),2*N+1,2*N+1),
			tmp,			2*N+1,2*N+1);		
	}
	
	if(par->verbose >0) {fprintf(stdout,"\n");}
	
	/* Libération de la mémoire */
	free(T11[0]); free(T11); free(T12[0]); free(T12);
	free(T21[0]); free(T21); free(T22[0]); free(T22);
	free(Z[0]); free(Z); free(tmp[0]); free(tmp);
	free(tmp2[0]); free(tmp2);
	free(tmp3[0]); free(tmp3);
	free(F_plus); free(F_moins); free(F_plus2); free(F_moins2);
	
	return 0;
}

/* Inutile ne marche pas ... problèmes numériques */
int md2D_comb_mat_S(complex **S11, complex **S12, complex **S21, complex **S22, 
		complex **S11_1, complex **S12_1, complex **S21_1, complex **S22_1, 
		complex **S11_2, complex **S12_2, complex **S21_2, complex **S22_2, int taille_matrice)
{
	int N = taille_matrice;
	complex **Q, **W, **Id, **tmp1, **tmp2, **tmp3;
	
	Q = allocate_CplxMatrix(N,N);
	W = allocate_CplxMatrix(N,N);
	Id = allocate_CplxMatrix(N,N);
	tmp1 = allocate_CplxMatrix(N,N);
	tmp2 = allocate_CplxMatrix(N,N);
	tmp3 = allocate_CplxMatrix(N,N);
	
	/* Q = inv [1 - S12_2 x S21_1 ]*/
	invM(Q,
		sub_M(tmp1,
			M_Id(Id,N),
			M_x_M(tmp2,
				S12_2,
				S21_1,	N,N),N,N),N);

	/* W = inv [1 - S21_1 x S12_2 ]*/
	invM(W,
		sub_M(tmp1,
			M_Id(Id,N),
			M_x_M(tmp2,
				S21_1,
				S12_2,	N,N),N,N),N);
	
	/* S11 = S11_1 x Q x S11_2 */
	M_x_M(S11,
		S11_1,
		M_x_M(tmp1,
			Q,
			S11_2,	N,N),N,N);	

	/* S22 = S22_2 x W x S22_1 */
	M_x_M(S22,
		S22_2,
		M_x_M(tmp1,
			W,
			S22_1,	N,N),N,N);	
	
	/* S12 = S12_1 + S11_1 x Q x S12_2 x S22_1 */
	add_M(S12,
		S12_1,
		M_x_M(tmp1,
			S11_1,
			M_x_M(tmp2,
				Q,
				M_x_M(tmp3,
					S12_2,
					S22_1,	N,N),N,N),N,N),N,N); 
	
	/* S21 = S21_2 + S22_2 x W x S21_1 x S11_2 */
	add_M(S21,
		S21_2,
		M_x_M(tmp1,
			S22_2,
			M_x_M(tmp2,
				W,
				M_x_M(tmp3,
					S21_1,
					S11_2,	N,N),N,N),N,N),N,N); 
	
	
	free(Q[0]);
	free(Q);
	free(W[0]);
	free(W);
	free(Id[0]);
	free(Id);
	free(tmp1[0]);
	free(tmp1);
	free(tmp2[0]);
	free(tmp2);
	free(tmp3[0]);
	free(tmp3);
	
	
	
	
	return 0;
}
