


/*-------------------------------------------------------------------------------------*/
/*!	\fn		int M_matrix(COMPLEX **M, double z, struct Param_struct *par)
 *
 *		\brief	M matrix calculation
 */
/*-------------------------------------------------------------------------------------*/
int M_matrix(COMPLEX **M, double z, struct Param_struct *par)
{
	int i,j;
	int vec_size = par->vec_size;

	COMPLEX *sigma_x = par->sigma_x;
	COMPLEX sigma_y0 = par->sigma_y0;
	COMPLEX **M11,**M13,**M14, **M21,**M23,**M24, **M31,**M32, **M41,**M42,**M43,**mM44;
	COMPLEX **invQzzQzx, **invQzzSigmax, **invQzzSigmay0, **Qtmp;
	COMPLEX **Qxx, **Qxz, **Qyy, **Qzx, **Qzz, **invQzz;
	Qxx = par->Qxx;
	Qxz = par->Qxz;
	Qyy = par->Qyy;
	Qzx = Qxz;
	Qzz = par->Qzz;
	invQzz = par->Qzz_1;
	
	M11 = allocate_CplxMatrix(vec_size,vec_size);
	M13 = allocate_CplxMatrix(vec_size,vec_size);
	M14 = allocate_CplxMatrix(vec_size,vec_size);
	M21 = allocate_CplxMatrix(vec_size,vec_size);
	M23 = allocate_CplxMatrix(vec_size,vec_size);
	M24 = allocate_CplxMatrix(vec_size,vec_size);
	M31 = allocate_CplxMatrix(vec_size,vec_size);
	M32 = allocate_CplxMatrix(vec_size,vec_size);
	M41 = allocate_CplxMatrix(vec_size,vec_size);
	M42 = allocate_CplxMatrix(vec_size,vec_size);
	M43 = allocate_CplxMatrix(vec_size,vec_size);
	mM44 = allocate_CplxMatrix(vec_size,vec_size);
	invQzzQzx = allocate_CplxMatrix(vec_size,vec_size);
	invQzzSigmax = allocate_CplxMatrix(vec_size,vec_size);
	invQzzSigmay0 = allocate_CplxMatrix(vec_size,vec_size);
	Qtmp = allocate_CplxMatrix(vec_size,vec_size);

	/* Toeplitz matrices calculations */
	mdC_QMatrix(z, Qxx, Qxz, Qyy, Qzz, invQzz, par);
	/* invQzzQzx = invQzz*Qzx */
	M_x_M(invQzzQzx, invQzz, Qzx, vec_size, vec_size);
	/* invQzzQzy = invQzz*Qzy */
	M_x_M(invQzzQzy, invQzz, Qzy, vec_size, vec_size);
	/* invQzzSigmax = invQzz*sigma_x */
	for(i=0;i<=vec_size-1;i++){
		for(j=0;j<=vec_size-1;j++){
			invQzzSigmax[i][j] = invQzz[i][j]*sigma_x[j];
		}
	}
	/* invQzzSigmay0 = invQzz*sigma_y0 */
	for(i=0;i<=vec_size-1;i++){
		for(j=0;j<=vec_size-1;j++){
			invQzzSigmay0[i][j] = invQzz[i][j]*sigma_y0;
		}
	}

	/* M11 = -sigma_x*invQzzQzx */
	for(i=0;i<=vec_size-1;i++){
		for(j=0;j<=vec_size-1;j++){
			M11[i][j] = -sigma_x[i]*invQzzQzx[i][j];
		}
	}
	/* M12 = 0 */

	/* M13 = sigma_x*invQzzSigmay0 */
	for(i=0;i<=vec_size-1;i++){
		for(j=0;j<=vec_size-1;j++){
			M13[i][j] = sigma_x[i]*invQzzSigmay0[i][j];
		}
	}
	/* M14 = Id - sigma_x*invQzzSigmax */
	for(i=0;i<=vec_size-1;i++){
		for(j=0;j<=vec_size-1;j++){
			M14[i][j] = -sigma_x[i]*invQzzSigmax[i][j];
		}
		M14[i][i] += 1;
	}

	/* M21 = -sigma_y0*invQzzQzx */
	for(i=0;i<=vec_size-1;i++){
		for(j=0;j<=vec_size-1;j++){
			M21[i][j] = -sigma_y0*invQzzQzx[i][j];
		}
	}
	/* M22 = 0 */

	/* M23 = -Id + sigma_y0*invQzzSigmay0 */
	for(i=0;i<=vec_size-1;i++){
		for(j=0;j<=vec_size-1;j++){
			M23[i][j] = sigma_y0*invQzzSigmay0[i][j];
		}
		M23[i][i] += -1;
	}
	/* M24 = -sigma_y0*invQzzSigmax */
	for(i=0;i<=vec_size-1;i++){
		for(j=0;j<=vec_size-1;j++){
			M24[i][j] = -sigma_y0*invQzzSigmax[i][j];
		}
	}

	/* M31 = -sigma_x*sigma_y0 */
	for(i=0;i<=vec_size-1;i++){
		for(j=0;j<=vec_size-1;j++){
			M31[i][j] = 0;
		}
		M31[i][i] -= sigma_x[i]*sigma_y0;
	}
	/* M32 = sigma_x^2 - Qyy */
	for(i=0;i<=vec_size-1;i++){
		for(j=0;j<=vec_size-1;j++){
			M32[i][j] = -Qyy[i][j];
		}
		M32[i][i] += sigma_x[i]*sigma_x[i];
	}

	/* M33 = 0 */

	/* M34 = 0 */

	/* M41 = -sigma_y0^2 + Qxx - Qxz*invQzzQzx */
	M_x_M(Qtmp,Qxz,invQzzQzx, vec_size, vec_size);
	for(i=0;i<=vec_size-1;i++){
		for(j=0;j<=vec_size-1;j++){
			M41[i][j] = Qxx[i][j] - Qtmp[i][j];
		}
		M41[i][i] -= sigma_y0*sigma_y0;
	}
	/* M42 = sigma_y*sigma_x */
	for(i=0;i<=vec_size-1;i++){
		for(j=0;j<=vec_size-1;j++){
			M42[i][j] = 0;
		}
		M42[i][i] += sigma_y0*sigma_x[i];
	}
	/* M43 = Qxz*invQzzSigmay0 */
	M_x_M(M43,Qxz,invQzzSigmay0, vec_size, vec_size);
	/* M44 = -Qxz*invQzzSigmax , mM44 = -M44 */
	M_x_M(mM44,Qxz,invQzzSigmax, vec_size, vec_size);

	/* putting it all together */
	for(i=0;i<=vec_size-1;i++){
		for(j=0;j<=vec_size-1;j++){
			M[i][j]            = I * M11[i][j];
			M[i][j+vec_size]   = 0;
			M[i][j+2*vec_size] = I * M13[i][j];
			M[i][j+3*vec_size] = I * M14[i][j];

			M[i+vec_size][j]            = I * M21[i][j];
			M[i+vec_size][j+vec_size]   = 0;
			M[i+vec_size][j+2*vec_size] = I * M23[i][j];
			M[i+vec_size][j+3*vec_size] = I * M24[i][j];

			M[i+2*vec_size][j]            =  I * M31[i][j];
			M[i+2*vec_size][j+vec_size]   =  I * M32[i][j];
			M[i+2*vec_size][j+2*vec_size] = 0;
			M[i+2*vec_size][j+3*vec_size] = 0;

			M[i+3*vec_size][j]            =  I * M41[i][j];
			M[i+3*vec_size][j+vec_size]   =  I * M42[i][j];
			M[i+3*vec_size][j+2*vec_size] =  I * M43[i][j];
			M[i+3*vec_size][j+3*vec_size] = -I * mM44[i][j];
		}
	}

#if 0
printf("\nz=%f\n",z);
printf("\nRe(M) :\n");SaveMatrix2file (M, 4*vec_size, 4*vec_size, "Re", "stdout");
printf("\nIm(M) :\n");SaveMatrix2file (M, 4*vec_size, 4*vec_size, "Im", "stdout");
printf("\nRe(M11) :\n");SaveMatrix2file (M11, vec_size, vec_size, "Re", "stdout");
printf("\nIm(M11) :\n");SaveMatrix2file (M11, vec_size, vec_size, "Im", "stdout");
printf("\nRe(M12) :\n");SaveMatrix2file (M12, vec_size, vec_size, "Re", "stdout");
printf("\nIm(M12) :\n");SaveMatrix2file (M12, vec_size, vec_size, "Im", "stdout");
printf("\nRe(M13) :\n");SaveMatrix2file (M13, vec_size, vec_size, "Re", "stdout");
printf("\nIm(M13) :\n");SaveMatrix2file (M13, vec_size, vec_size, "Im", "stdout");
printf("\nRe(M14) :\n");SaveMatrix2file (M14, vec_size, vec_size, "Re", "stdout");
printf("\nIm(M14) :\n");SaveMatrix2file (M14, vec_size, vec_size, "Im", "stdout");
printf("\nRe(M21) :\n");SaveMatrix2file (M21, vec_size, vec_size, "Re", "stdout");
printf("\nIm(M21) :\n");SaveMatrix2file (M21, vec_size, vec_size, "Im", "stdout");
printf("\nRe(M22) :\n");SaveMatrix2file (M22, vec_size, vec_size, "Re", "stdout");
printf("\nIm(M22) :\n");SaveMatrix2file (M22, vec_size, vec_size, "Im", "stdout");
printf("\nRe(M23) :\n");SaveMatrix2file (M23, vec_size, vec_size, "Re", "stdout");
printf("\nIm(M23) :\n");SaveMatrix2file (M23, vec_size, vec_size, "Im", "stdout");
printf("\nRe(M24) :\n");SaveMatrix2file (M24, vec_size, vec_size, "Re", "stdout");
printf("\nIm(M24) :\n");SaveMatrix2file (M24, vec_size, vec_size, "Im", "stdout");
printf("\nRe(M31) :\n");SaveMatrix2file (M31, vec_size, vec_size, "Re", "stdout");
printf("\nIm(M31) :\n");SaveMatrix2file (M31, vec_size, vec_size, "Im", "stdout");
printf("\nRe(M32) :\n");SaveMatrix2file (M32, vec_size, vec_size, "Re", "stdout");
printf("\nIm(M32) :\n");SaveMatrix2file (M32, vec_size, vec_size, "Im", "stdout");
printf("\nRe(mM33) :\n");SaveMatrix2file (mM33, vec_size, vec_size, "Re", "stdout");
printf("\nIm(mM33) :\n");SaveMatrix2file (mM33, vec_size, vec_size, "Im", "stdout");
printf("\nRe(M34) :\n");SaveMatrix2file (M34, vec_size, vec_size, "Re", "stdout");
printf("\nIm(M34) :\n");SaveMatrix2file (M34, vec_size, vec_size, "Im", "stdout");
printf("\nRe(M41) :\n");SaveMatrix2file (M41, vec_size, vec_size, "Re", "stdout");
printf("\nIm(M41) :\n");SaveMatrix2file (M41, vec_size, vec_size, "Im", "stdout");
printf("\nRe(M42) :\n");SaveMatrix2file (M42, vec_size, vec_size, "Re", "stdout");
printf("\nIm(M42) :\n");SaveMatrix2file (M42, vec_size, vec_size, "Im", "stdout");
printf("\nRe(M43) :\n");SaveMatrix2file (M43, vec_size, vec_size, "Re", "stdout");
printf("\nIm(M43) :\n");SaveMatrix2file (M43, vec_size, vec_size, "Im", "stdout");
printf("\nRe(mM44) :\n");SaveMatrix2file (mM44, vec_size, vec_size, "Re", "stdout");
printf("\nIm(mM44) :\n");SaveMatrix2file (mM44, vec_size, vec_size, "Im", "stdout");
#endif

	free(M11[0]);
	free(M13[0]);
	free(M14[0]);
	free(M21[0]);
	free(M23[0]);
	free(M24[0]);
	free(M31[0]);
	free(M32[0]);
	free(M41[0]);
	free(M42[0]);
	free(M43[0]);
	free(mM44[0]);
	free(invQzzQzx[0]);
	free(invQzzSigmax[0]);
	free(invQzzSigmay0[0]);
	free(Qtmp[0]);
	free(M11);
	free(M13);
	free(M14);
	free(M21);
	free(M23);
	free(M24);
	free(M31);
	free(M32);
	free(M41);
	free(M42);
	free(M43);
	free(mM44);
	free(invQzzQzx);
	free(invQzzSigmax);
	free(invQzzSigmay0);
	free(Qtmp);
	
	return 0;
}

/*-------------------------------------------------------------------------------------*/
/*!	\fn		int QMatrix(double z, COMPLEX **Qxx, COMPLEX **Qxy, COMPLEX **Qxz, COMPLEX **Qyy, COMPLEX **Qyz, COMPLEX **Qzz, COMPLEX **Qzz_1, struct Param_struct *par)
 *
 *		\brief	Q matrix calculation
 */
/*-------------------------------------------------------------------------------------*/
int mdC_QMatrix(double z, COMPLEX **Qxx, COMPLEX **Qxz, COMPLEX **Qyy, COMPLEX **Qzz, COMPLEX **Qzz_1, struct Param_struct *par)
{
	int vec_size=par->vec_size, N_x=par->N_x, N=par->N;
	COMPLEX **M_buff1, **toepk2, **invtoepk2, **minusDelta, *buffer_N_xNpry;

	M_buff1 = allocate_CplxMatrix(vec_size,vec_size);
	toepk2 = allocate_CplxMatrix(vec_size,vec_size);
	invtoepk2 = allocate_CplxMatrix(vec_size,vec_size);
	minusDelta = allocate_CplxMatrix(vec_size,vec_size);

	/* Calculating normal vectors Fourier Toeplitz matrices */
	mdC_normalCoef(z, par);

	buffer_N_x = (COMPLEX *) malloc(sizeof(COMPLEX)*N_x);


	/* toepk2 = [[k2]] */
	toeplitz_2D(toepk2, N, (*par->k_2)(par, buffer_N_x, z), N_x, Npry);
	toeplitz_1D(toepk2, N, COMPLEX *f_x, int N_x)

	/* invtoepk2 = inv([[1/k2]]) */
	toeplitz_2D(M_buff1, Nx, Ny, (*par->invk_2)(par, buffer_Nx, z), N_x, Npry);
	invM(invtoepk2, M_buff1, vec_size);




	/* minusDelta = inv([[1/k2]]) - [[k2]] */
	sub_M(minusDelta, invtoepk2, toepk2, vec_size, vec_size);

	/* Qxx = toepk2 -DeltaNxx */
	add_M(Qxx, toepk2, M_x_M(M_buff1, minusDelta, par->Nxx, vec_size, vec_size), vec_size, vec_size);

	/* Qxz = -DeltaNxz */
	M_x_M(Qxz,minusDelta,par->Nxz,vec_size,vec_size);

	/* Qyy = toepk2 */
	M_equal(Qyy, toepk2, vec_size, vec_size);

	/* Qzz = toepk2 -DeltaNzz */
	add_M(Qzz, toepk2, M_x_M(M_buff1,minusDelta,par->Nzz,vec_size,vec_size), vec_size, vec_size);

	/* Qzz_1 = inv(Qzz) */
	invM(Qzz_1, Qzz, vec_size);
	
#if 0
printf("\n Rek2\n");SaveCplxTab2file ((*par->k_2)(par, buffer_N_xNpry, z), Npry*N_x, "Re", "stdout", " ", N_x, "\n");
SaveMatrix2file (toepk2, vec_size, vec_size, "Re", "stdout");
printf("\n toepk2\n");SaveMatrix2file (toepk2, vec_size, vec_size, "Re", "stdout");
printf("\n invtoepk2\n");SaveMatrix2file (invtoepk2, vec_size, vec_size, "Re", "stdout");
printf("\n moinsDelta\n");SaveMatrix2file (moinsDelta, vec_size, vec_size, "Re", "stdout");
printf("\n Qxx\n");SaveMatrix2file (Qxx, vec_size, vec_size, "Re", "stdout");
printf("\n Qxy\n");SaveMatrix2file (Qxy, vec_size, vec_size, "Re", "stdout");
printf("\n Qxz\n");SaveMatrix2file (Qyz, vec_size, vec_size, "Re", "stdout");
printf("\n Qyy\n");SaveMatrix2file (Qyy, vec_size, vec_size, "Re", "stdout");
printf("\n Qyz\n");SaveMatrix2file (Qyz, vec_size, vec_size, "Re", "stdout");
printf("\n Qzz\n");SaveMatrix2file (Qzz, vec_size, vec_size, "Re", "stdout");
#endif
	free(buffer_Nx);
	free(M_buff1[0]);
	free(M_buff1);
	free(toepk2[0]);
	free(toepk2);
	free(invtoepk2[0]);
	free(invtoepk2);
	free(moinsDelta[0]);
	free(moinsDelta);
	
	return 0;
}

/*-------------------------------------------------------------------------------------*/
/*!	\fn		int toeplitz_1D (COMPLEX **M_toep, int N, COMPLEX *f_x, int N_x)
 *
 *		\brief	Returns the Toeplitz matrix of the (-N,+N) Fourier coefficients of function f_x

 */
/*-------------------------------------------------------------------------------------*/
int toeplitz_1D (COMPLEX **M_toep, int N, COMPLEX *f_x, int N_x)
{

	int i,j;
	int N_tf = 2*N;
	int vec_size = 2*N+1;
	fftw_plan plan_TF;
	
	COMPLEX *tf_long,*TF_1D;
	tf_long = (COMPLEX *) malloc(sizeof(COMPLEX)*N_x);
	TF_1D =   (COMPLEX *) malloc(sizeof(COMPLEX)*(4*N+1));
	
	/* Computing the FFT */
	plan_TF = fftw_plan_dft_1d(N_x, (fftw_complex *)f_x, (fftw_complex *)tf_long, FFTW_FORWARD, FFTW_ESTIMATE);	
	fftw_execute(plan_TF);

	/* Keeping only the components between -N_tf & +N_tf and normalizing by 1/N_x */
	double coefnorm = 1.0/N_x;
	for (i=0;i<=N_tf-1;i++){
		TF_1D   [i]      = tf_long [N_x-N_tf+i] * coefnorm;
		TF_1D   [i+N_tf] = tf_long [i]          * coefnorm;
	}
	TF_1D      [2*N_tf] = tf_long [N_tf]       * coefnorm;

	/* Building the Toeplitz matrix */
	for(i=0;i<=vec_size-1;i++){
		for(j=0;j<=vec_size-1;j++){
			M_toep[i][j] = TF_1D[vec_size-1+i-j];
		}
	}
	/* Freeing memory */
	fftw_destroy_plan(plan_TF);
	free(tf_long);	free(TF_1D);

	return 0;
}
