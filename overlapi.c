#include <stdio.h>

#include <math.h>
#include <stdlib.h>
#include "/mnt/c/ruby_home/physics.h" 

#define size 256 

#define m (126.90447*amu)
#define	mu (m*m/(m+m)) //ヨウ素の換算質量	

#define	k (pow((2*(M_PI)*c_spd*150*100),2)*mu) //力の定数 SIに変換してる
#define	a (2*(M_PI)*sqrt(k*mu)/plank) //アルファのこと
#define	d (pow((a/(M_PI)),0.25))


double dsyevx_(
		char *JOBZ;
		char *RANGE;
		char *UPLO;
		int *N;
		double *A;
		int *LDA;
		double *VL;
		double *VU;
		int *IL;
		int *IU;
		double *ABSTOL;
		int *M;
		double *W;
		double *Z;
		int *LDZ;
		double *WORK;
		int *LWORK;
		int *IWORK;
		int *IFAIL;
		int *INFO;

);

int main(void){
	int n,i,j;
	double dx,x0,xf,x,kc,S,s;
	double H[size][size];
	FILE *fp,*fp2;
	
	
	if((fp=fopen("diagonal.dat","w"))==NULL){
		fprintf(stderr,"ERROR\n");
		return EXIT_FAILURE;
	}	
	if((fp2=fopen("vector_0.dat","w"))==NULL){
		fprintf(stderr,"ERROR\n");
		return EXIT_FAILURE;
	}	
	n=size;
	char JOBZ='V';
	char RANGE='A';
	char UPLO='U';
	int N=size;
	double A[n][n];
	int M=size;
	int INFO;
	int LDA=n;
	int LDZ=n;
	int LWORK=8*n;
	double ABSTOL=1.0e-40;
	int IFAIL[n];
	int IWORK[5*n];
	double W[n];
	double WORK[LWORK];
	double Z[LDA][n];
	double VL;
	double VU;
	int IL;
	int IU;


	x0=-1.0e-10;
	xf=1.0e-10;
	dx=(xf-x0)/n;
	n=size;
	

//hamiltonian
	kc=plank*plank/(8*(M_PI)*(M_PI)*mu*dx*dx);
	for(i=0; i<=n-1; i++){
		x=x0+dx*i;	
		for(j=0; j<=n-1; j++){
			
			if(i==j){
				H[i][j]=kc*(M_PI)*(M_PI)/(3)+k*x*x/2;
			}	
			else{
				H[i][j]=kc*2.0/((i-j)*(i-j))*pow(-1,(double)(i-j));
			}
		}	
	}


	for(i=0;i<size;i++){
		for(j=0; j<size; j++){
			A[i][j]=H[j][i];
		}
	}
	
dsyevx_(&JOBZ, &RANGE, &UPLO, &N, A, &LDA, &VL, &VU, &IL, &IU, &ABSTOL, &M, W, Z, &LDZ, WORK, &LWORK, IWORK, IFAIL, &INFO);
/*	printf("%1.15e\n", H[1][1]);
	printf("%1.15e\n", H[1][0]);
	printf("%1.15e\n", k);
	printf("%1.15e\n", kc);
*/
	printf("%1.15e\n",a);
	
	double ps0[n];
	double ps1[n];
	double ps2[n];
	double ps3[n];
	double ps4[n];


	for(i=0;i<=n-1;i++){
		x=x0+dx*i;	

		ps0[i]=(pow(a/M_PI,0.25)*exp(-a*x*x/2))*sqrt(dx);
		ps1[i]=(pow(4*a*a*a/M_PI,0.25)*x*exp(-a*x*x/2))*sqrt(dx);
		ps2[i]=(pow(a/(4*M_PI),0.25)*(2*a*x*x-1)*exp(-a*x*x/2))*sqrt(dx);
		ps3[i]=(pow(a*a*a/(9*M_PI),0.25)*(2*a*x*x*x-3*x)*exp(-a*x*x/2))*sqrt(dx);
		ps4[i]=(sqrt(1.0/24)*pow(a/(M_PI),0.25)*(4*a*a*x*x*x*x-12*a*x*x+3)*exp(-a*x*x/2))*sqrt(dx);
	}

	for(i=0; i<size; i++){
		
		x=x0+dx*i;	
			fprintf(fp,"%1.15e\n",W[i]/(100*plank*c_spd));		
			fprintf(fp2,"%1.15e %1.15e %1.15e\n",x,Z[0][i],ps0[i]);		

	}
	fclose(fp);
	fclose(fp2);

	S=0;
	for(i=0;i<=n-1;i++){
		x=x0+dx*i;	
		s=ps0[i]*Z[0][i];
	//	s=Z[0][i]*Z[0][i];
		S=S+s;
	
	}
	printf("S0=%1.15e\n",S );

	S=0;
	for(i=0;i<=n-1;i++){
		x=x0+dx*i;	
		s=ps1[i]*Z[1][i];
		S=S+s;
	}	
	printf("S1=%1.15e\n",S );

	S=0;
	for(i=0;i<=n-1;i++){
		x=x0+dx*i;	
		s=ps2[i]*Z[2][i];
		S=S+s;
	}	
	printf("S2=%1.15e\n",S );

	S=0;
	for(i=0;i<=n-1;i++){
		x=x0+dx*i;	
		s=ps3[i]*Z[3][i];
		S=S+s;
	}	
	printf("S3=%1.15e\n",S );

	S=0;
	for(i=0;i<=n-1;i++){
		x=x0+dx*i;	
		s=ps4[i]*Z[4][i];
		S=S+s;
	}	
	printf("S4=%1.15e\n",S );

return 0;

}








