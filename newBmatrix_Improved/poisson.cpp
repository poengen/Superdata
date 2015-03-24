#include <stddef.h>
#include <stdlib.h>
#include <stdio.h>
#include <memory.h>
#include <math.h>
#include <cmath>
#include <iostream>
#include <mpi.h>
#include <time.h>
typedef double Real;
using namespace std;

/* function prototypes */
Real *createRealArray (int n);
Real **createReal2DArray (int m, int n);
void transpose (Real **bt, Real **b, int n, int m);
extern "C"{
void fst_(Real *v, int *n, Real *w, int *nn);
void fstinv_(Real *v, int *n, Real *w, int *nn);
}

//print functions
void printArray(Real **A, int n1, int n2);

// andre funksjoner
int *createIntArray (int n);


int main(int argc, char **argv )
{
	//declaring variables
	//clock_t start, end; start = clock();
  	Real *diag, **b, **bt, *z, **new_b , **analyticSol, *bPointer, *ErrorVector;
  	Real pi, h, umax;
  	int i, j, n, m, nn;
  	
  	int *sendcounts, *displs;
        div_t divresult; 
        int local_col, displace;

	//initializing MPI
	int rank, size, k , l, current, prev;
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);


  /* the total number of grid points in each spatial direction is (n+1) */
  /* the total number of degrees-of-freedom in each spatial direction is (n-1) */
  /* this version requires n to be a power of 2 */

 	if( argc < 2 ) { printf("need a problem size\n"); return 1; }

  	n  = atoi(argv[1]);
  	m  = n-1;
  	nn = 4*n;
  	
  	h    = 1./(Real)n;
  	pi   = 4.*atan(1.);
  	
  	divresult = div(m,size);
  	local_col = divresult.quot;
  	if (rank < divresult.rem){ local_col++;}
	
  	diag = createRealArray (m);
  	b    = createReal2DArray (local_col,m);
	analyticSol = createReal2DArray (local_col,m);	
  	bt   = createReal2DArray (m,local_col);
  	z    = createRealArray (nn);
  	
  	
	// Lager en int vektor med antall elementer som skal sendes til hver prosessor. 
	sendcounts =  createIntArray(size);
	for (i = 0; i < size; i++){
		if (i < divresult.rem){ sendcounts[i] = local_col*(divresult.quot+1); }
		else { sendcounts[i]=local_col*divresult.quot; }
	}
	
	

	// lager en displacement vektor som
	displs = createIntArray(size);
	displs[0]=0;
	for (i = 1 ; i<size ; i++){ displs[i]=displs[i-1]+sendcounts[i-1]; }

	
	// Finner den globale indexen til den foerste kolonna
	displace = 0; 
	for (i = 0 ; i < (rank); i++){
		if (i<divresult.rem){ displace+=(divresult.quot+1);}
		else { displace+=divresult.quot;}
	}	


  	for (i=0; i < m; i++) { diag[i] = 2.*(1.-cos((i+1)*pi/(Real)n)); }
  	
  	
  	// Lager b-matrisa og matrisa som inneholder eksakt loesning paa problemet.
  	for (j=0; j < local_col; j++) {
    		for (i=0; i < m; i++) {
			analyticSol[j][i]= sin(pi*h*(i+1))*sin(2*pi*h*(1+displace+j));
      			b[j][i] =h*h*5*pi*pi*sin(pi*h*(i+1))*sin(2*pi*h*(1+displace+j));
    		}
  	}
  	
  	for (j=0; j < local_col; j++) { fst_(b[j], &n, z, &nn); }
  	
	
  	//// FORSTE TRANSPONERING ////
  	// Transponerer matrisen loakt
  	transpose(bt,b,m,local_col);
  	//Sender info
	MPI_Alltoallv(bt[0],sendcounts,displs,MPI_DOUBLE, b[0], sendcounts,displs,MPI_DOUBLE, MPI_COMM_WORLD );
	// Prosedyren som legger motatt info i riktig rekkefoelge
	new_b    = createReal2DArray (local_col,m);
	bPointer = b[0];
	prev = 0 ; current = 0; 
	for (i = 0; i <size; i++){
		 l = -1 ; 
		prev += current; 
		if (i < divresult.rem ){ current = divresult.quot+1;}
		else { current = divresult.quot; }
		for (j = 0 ; j < sendcounts[i]; j++){
			if (j%current ==0 ){ l++; k=prev; }
			new_b[l][k]=bPointer[0];
			bPointer++; k++;
		}
	}
	//// SLUTT FORSTE TRANSPONERING ////

	for (i = 0 ; i < local_col ; i++) {fstinv_(new_b[i],&n,z,&nn); }
	
	
	for ( j = 0; j <local_col; j++){
		for (i = 0 ; i < m ; i++){
			new_b[j][i]=new_b[j][i]/(diag[i]+diag[j+displace]);
		}
	}
	
	for (i = 0 ; i < local_col ; i++){ fst_(new_b[i],&n,z,&nn); } 
	
	
	//// ANDRE TRANSPONERING /////
	transpose(bt,new_b,m,local_col);
	MPI_Alltoallv(bt[0],sendcounts,displs,MPI_DOUBLE, b[0], sendcounts,displs,MPI_DOUBLE, MPI_COMM_WORLD );
	bPointer = b[0]; k = 0 ; l = 0; current = 0 ; prev = 0;
	for (i = 0; i <size; i++){
		 l = -1 ; 
		prev += current; 
		if (i < divresult.rem ){ current = divresult.quot+1;}
		else { current = divresult.quot; }
		for (j = 0 ; j < sendcounts[i]; j++){
			if (j%current ==0 ){ l++; k=prev; }
			new_b[l][k]=bPointer[0];
			bPointer++; k++;
		}
	}
	//// SLUTT ANDRE TRANSPONERING /////

	for (j = 0 ; j < local_col; j++ ){ fstinv_(new_b[j],&n,z,&nn); }
	
	
	//// FINNER ERROR /////
        umax = 0;
  	for (j=0; j < local_col; j++) {
    		for (i=0; i < m; i++) {
			h = abs(new_b[j][i]-analyticSol[j][i]);	
      			if (h > umax) {umax = h;}
    		}
  	}

	if (rank==0){ ErrorVector = createRealArray(size); }
	else { ErrorVector=NULL; }
        // Sender max error fra hver prosessor til prosessor 0
	MPI_Gather(&umax,1,MPI_DOUBLE,ErrorVector,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
        
        if (rank==0){ 
                for (i = 1 ; i < size ; i ++){
                if (ErrorVector[i]>umax) umax=ErrorVector[i];
                }
                cout<<" umax= "<<umax<<endl;
        }

        
  	
  	
  	
  	
  	
  	MPI_Finalize();
  	return 0;
}

////////////////////////////////////////////////////////////
void transpose (Real **bt, Real **b, int n, int m)
{
  int i, j;
  for (j=0; j < n; j++) {
    for (i=0; i < m; i++) {
      bt[j][i] = b[i][j];
    }
  }
}

Real *createRealArray (int n)
{
  Real *a;
  int i;
  a = (Real *)malloc(n*sizeof(Real));
  for (i=0; i < n; i++) {
    a[i] = 0.0;
  }
  return (a);
}

int *createIntArray (int n){
  int *a;
  int i;
  a = (int *)malloc(n*sizeof(int));
  for (i=0; i < n; i++) {
    a[i] = 0.0;
  }
  return (a);
	
}

Real **createReal2DArray (int n1, int n2)
{
  int i, n;
  Real **a;
  a    = (Real **)malloc(n1   *sizeof(Real *));
  a[0] = (Real  *)malloc(n1*n2*sizeof(Real));
  for (i=1; i < n1; i++) {
    a[i] = a[i-1] + n2;
  }
  n = n1*n2;
  memset(a[0],0,n*sizeof(Real));
  return (a);
}

void printArray(Real **A, int n1, int n2){
	cout<<endl;
	for (int i = 0 ; i < n1; i++){
		for (int j = 0; j < n2; j++){
			cout << A[i][j]<<" ";
		}
		cout<<endl;
	}
	cout<<endl;
}




















