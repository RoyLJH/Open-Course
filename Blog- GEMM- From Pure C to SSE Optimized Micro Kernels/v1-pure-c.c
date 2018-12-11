#include <stdio.h>
#include <memory.h>
/*
    Pay Attention:
	All the original matrix are saved column by column,
	we pack A into MR * KC * (MC/MR) column-wise matrix
	and pack B into KC * NR * (NC/NR) row-wise matrix
*/
#define M 16
#define K 15
#define N 16

#define MC 8
#define NC 12
#define KC 11
#define MR 4
#define NR 4

double  A[M*K];
double  B[K*N];
double  C[M*N];

double _A[MC*KC];
double _B[KC*NC];
double _C[MR*NR];

void
initMatrix(int m, int n, double *X, int ldX, int counter)
{
	//your code here
	int i;
	for(i=0; i<n*m; i++){
	    X[i] = counter++;
	}
}

void
printMatrix(int m, int n, const double *X, int ldX)
{
	//your code here
	int r,c;
	double v;	
        for(r=0; r<m; r++){
            for(c=0; c<n; c++){
	        v = X[c*ldX+r];
	        //if(v<10) printf("   ");
    	        //else if(v<100) printf("  ");
	        //else printf(" ");
	        printf("  %-8.2lf",v);
	    }
	    printf("\n");
	}
	
}

void
pack_MRxk(int k, const double *A, int incRowA, int incColA, double *buffer)
{
	//your code here
	int i, j;
	for(i=0; i<k; i++){
	    for(j=0; j<MR; j++){
		buffer[j] = A[j*incRowA];
	    }
	    buffer += MR;
	    A += incColA;
	}
	
}

void
pack_A(int mc, int kc, const double *A, int incRowA, int incColA, double *buffer)
{
	//your code here
	memset(buffer, 0, mc*kc*sizeof(double));
	int panels = mc/MR;
	int rows_remain = mc%MR;
	int p, c;
	for(p=0; p<panels; p++){
	    pack_MRxk(kc, A, incRowA, incColA, buffer);
	    A += MR*incRowA;
	    buffer += MR*kc;
	}
	if(rows_remain > 0){
	    for(c=0; c<kc; c++){
	    	for(p=0; p<rows_remain; p++){
		    buffer[p] = A[p*incRowA];
		}
		for(; p<MR; p++){
		    buffer[p] = 0.0;
		}
		buffer += MR;
		A += incColA;
	    }
	}
}

void
pack_kxNR(int k, const double *B,int incRowB, int incColB, double *buffer){
	int i,j;
	for(i=0; i<k; i++){
	    for(j=0; j<NR; j++){
		buffer[j] = B[j*incColB];
	    }
	    buffer += NR;
	    B += incRowB;
	}
}

void
pack_B(int kc, int nc, const double*B, int incRowB, int incColB, double *buffer){
	memset(buffer, 0, nc*kc*sizeof(double));
	int panels = nc/NR;
	int cols_remain = nc%NR;
	int p, r;
	for(p=0; p<panels; p++){
	    pack_kxNR(kc, B, incRowB, incColB, buffer);
	    buffer += kc*NR;
	    B += NR*incColB;
	}
	if(cols_remain>0){
	    for(r=0; r<kc; r++){
		for(p=0; p<cols_remain; p++){
		    buffer[p] = B[p*incColB];
		}
		for(; p<NR; p++){
		    buffer[p] = 0.0;
		}
		buffer += NR;
		B += incRowB;
	    }
	}
}

// Multiplying one panel from _A and one panel from _B

void
dgemm_micro_kernel(int kc, double alpha, const double*A, const double*B, double beta,
		double*C, int incRowC, int incColC)
{
	// your code here, implementing C <- beta*C + alpha*A*B, where C[MR*NR]
	double AB[MR*NR];
	int m,n,k;

	// AB = A*B
	memset(AB, 0, MR*NR*sizeof(double));
	
	for(k=0; k<kc; k++){
	    for(m=0; m<MR; m++){
		for(n=0; n<NR; n++){
		    AB[n*MR+m] += A[m] * B[n];
		}
	    }
	    A += MR;
	    B += NR;
	}
	
	// C <- beta*C
	if(beta == 0.0){
	    //Pay Attention: Do Not Use memset to Initialize C because this MR*NR area is not memory-continuous!!
	    //memset(C, 0, MR*NR*sizeof(double));
	    for(m=0; m<MR; m++){
		for(n=0; n<NR; n++){
		    C[n*incColC+m*incRowC] = 0.0;
		}
	    }
	}else if(beta != 1.0){
	    for(m=0; m<MR; m++){
		for(n=0; n<NR; n++){
		    C[n*incColC+m*incRowC] *= beta;
		}
	    }
	}

	// C <- C + alpha*AB

	if(alpha == 0.0){
	}else if(alpha == 1.0){
	    for(m=0; m<MR; m++){
		for(n=0; n<NR; n++){
		    C[n*incColC+m*incRowC] += AB[n*MR+m];
		}
	    }
	}else{
	    for(m=0; m<MR; m++){
		for(n=0; n<NR; n++){
		    C[n*incColC+m*incRowC] += alpha*AB[n*MR+m];
		}
	    }
	}
}

// Multiplication of blocks of A and B (Composed of several panels), which were previously packed to buffers _A and _B

void dgemm_macro_kernel(int mc,int nc,int kc,double alpha,double beta,
			double*C,int incRowC,int incColC)
{
	int mp = (mc+MR-1)/MR;
	int np = (nc+NR-1)/NR;
	int rows_remain = mc%MR;
	int cols_remain = nc%NR;
	int m, n, r, c, mr, nr;
	for(m=0; m<mp; m++){
	    for(n=0; n<np; n++){
		if( (m!=mp-1 || rows_remain == 0) && (n!=np-1 || cols_remain == 0) ){	// If both panels are MR and NR full

		    dgemm_micro_kernel(kc, alpha, &_A[m*MR*kc], &_B[n*NR*kc], beta, &C[m*MR*incRowC+n*NR*incColC], incRowC, incColC);

		}else{	// If one of the panels is not full in MR or NR, namely have 0s

		    dgemm_micro_kernel(kc, 0.0, &_A[m*MR*kc], &_B[n*NR*kc], 1.0, _C, 1, MR);	// Can we reduce some calculations by neglecting 0 * x ?
		    
		    mr = (m!=mp-1 || rows_remain==0) ? MR : rows_remain;
		    nr = (n!=np-1 || cols_remain==0) ? NR : cols_remain;
		    for(r=m*MR; r<m*MR+mr; r++){
			for(c=n*NR; c<n*NR+nr; c++){
			    if(beta==0.0){
				C[r*incRowC+c*incColC] = 0.0;
			    }else if(beta!=1.0){
				C[r*incRowC+c*incColC] *= beta;
			    }
			    if(alpha==1.0){
				C[r*incRowC+c*incColC] += _C[(c-n*NR)*MR+(r-m*MR)];
			    }else if(alpha!=0.0){
				C[r*incRowC+c*incColC] += _C[(c-n*NR)*MR+(r-m*MR)]*alpha;
			    }
			}
		    }
		}
	    }
	}
}

void
dgemm_nn(int m, int n, int k, double alpha, const double *A, int incRowA, int incColA,
	const double *B, int incRowB, int incColB, double beta, double *C, int incRowC, int incColC)
{
	int mb = (m+MC-1)/MC;
	int nb = (n+NC-1)/NC;
	int kb = (k+KC-1)/KC;
	
	int m_r = m%MC;
	int n_r = n%NC;
	int k_r = k%KC;

	int mc, nc, kc;
	int i, j, l;

	double _beta;

	if(m==0 || k==0 || n==0 || alpha==0.0){		//no need for A*B, just C <- beta*C
	    if(beta != 1.0){
		for(i=0; i<M; i++){
		    for(j=0; j<N; j++){
			C[i*incRowC+j*incColC] *= beta;
		    }
		}
	    } 
	    return;
	}

	for(l=0; l<kb; l++){
	    kc = (l!=kb-1 || k_r==0)? KC : k_r;
	    _beta = (l==0)? beta : 1.0;
	    for(j=0; j<nb; j++){
		nc = (j!=nb-1 || n_r==0)? NC : n_r;
		pack_B(kc, nc, &B[l*KC*incRowB+j*NC*incColB], incRowB, incColB, _B);

		for(i=0; i<mb; i++){
		    mc = (i!=mb-1 || m_r==0) ? MC : m_r;
		    pack_A(mc, kc, &A[i*MC*incRowA+l*KC*incColA], incRowA, incColA, _A);

		    dgemm_macro_kernel(mc, nc, kc, alpha, _beta, &C[i*MC*incRowC+j*NC*incColC], incRowC, incColC);
		}
	    }
	}
}


int
main()
{
	//int i, j, mc, kc, nc;
	
	initMatrix(M, K, A, M, 1);
	initMatrix(K, N, B, K, M*K+1);

	printf("A = \n");
	printMatrix(M, K, A, M);

	printf("B = \n");
	printMatrix(K, N, B, K);

	//Checked the answer is right
	dgemm_nn(M, N, K, 1.0, A, 1, M, B, 1, K, 0.0, C, 1, M);

	printf("C = \n");
	printMatrix(M, N, C, M);

	return 0;
}

