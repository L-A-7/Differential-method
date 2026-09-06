% my first octave script

L=16000;
N=100;
n_super = 1.0 + 0.0i;
n_sub   = 1.5 + 0.0i;
lambda = 632.8
k_super = 2*pi*n_super/lambda
k_sub   = 2*pi*n_sub/lambda

delta_sigma = 2*pi/L;


% Loading reflected field
load "chp_out.txt" F
Ar=zeros(2*N+1,1);
for n=-N:N
	Ar(n+N+1) = F(1,n+N+1) + i*F(2,n+N+1);
endfor

% Calculating the Inverse Reflexion Matrix
inv_r = zeros(2*N+1,2*N+1);
for n=-N:N
	sigma = n*delta_sigma;
	% TE case
	b_super = sqrt(k_super^2-sigma^2);
	b_sub   = sqrt(k_sub^2  -sigma^2);
	r(n+N+1,n+N+1) = (b_super - b_sub) / (b_super + b_sub);
	inv_r(n+N+1,n+N+1) = (b_super + b_sub) / (b_super - b_sub);
endfor

% Calculating the inverse incident field
Ai = inv_r * Ar;

% Writing to a file
fp = fopen ("fileout.dat", "w", "ieee-le");
fprintf(fp,"Re_Ai = ");
for n=-N:N
	fprintf(fp,"%1.8e ",real(Ai(n+N+1)))
endfor
fprintf(fp,"\nIm_Ai = ");
for n=-N:N
	fprintf(fp,"%1.8e ",imag(Ai(n+N+1)))
endfor
fprintf(fp,"\n");
fclose(fp);

