#include <stddef.h>
#include <stdlib.h>
#include <stdio.h>
#include <memory.h>
#include <math.h>
#include <iostream>
#include <mpi.h>
#include <time.h>
typedef double Real;
using namespace std;

/* function prototypes */
Real *createRealArray (int n);
Real **createReal2DArray (int m, int n);
void transpose (Real **bt, Real **b, int m);
extern "C"{
void fst_(Real *v, int *n, Real *w, int *nn);
void fstinv_(Real *v, int *n, Real *w, int *nn);
}

//print functions
void print2DArray(int m, int local_quot, Real **b, int rank);
void print1DArray(int size, int array[], int rank);
void printArray(Real **A, int n1, int n2);

// andre funksjoner
int *createIntArray (int n);

int main(int argc, char **argv)
{
	//declaring variables
	//clock_t start, end; start = clock();
  	Real *diag, **b, **bt, *z, **temp_b;
  	Real pi, h, umax;
  	int i, j, n, m, nn;


	//initializing MPI
	int rank, size;
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);


  /* the total number of grid points in each spatial direction is (n+1) */
  /* the total number of degrees-of-freedom in each spatial direction is (n-1) */
  /* this version requires n to be a power of 2 */

 	if( argc < 2 ) {
    		printf("need a problem size\n");
    		return 1;
  	}

  	n  = atoi(argv[1]);
  	m  = n-1;
  	nn = 4*n;

	h = 1./(Real)n;
  	pi = 4.*atan(1.);

//up to date here
  	
  	div_t divresult = div(m,size);
  	
	int local_quot = divresult.quot; //divresult.quot is the same on all threads
  	if (rank < divresult.rem){
		local_quot++; //local_quot iterated by 1 on all divresult.rem first threads
	}

	//print rank and local_quot
	//cout << "Rank = " << rank << ", local_quot = " << local_quot << endl;

  	diag = createRealArray (m);
  	b    = createReal2DArray (local_quot,m);
  	
  	bt   = createReal2DArray (m,local_quot);
  	z    = createRealArray (nn);


  	for (i=0; i < m; i++) {
    		diag[i] = 2.*(1.-cos((i+1)*pi/(Real)n));
  	}

  	for (j=0; j < local_quot; j++) {
    		for (i=0; i < m; i++) {
      			b[j][i] = h*h;
    		}
  	}

	for (j=0; j < local_quot; j++) {
    		fst_(b[j], &n, z, &nn);
  	}



//denne funksjonen lager et array hvor array[i] er antall kolonner av b på prosess i
//den er lik på alle prosesser og kan brukes i reshape()-funksjonen
	int size_array[size];
	int temp = divresult.rem;
	for (i=0; i<size; i++){
		if (temp > 0){size_array[i] = divresult.quot + 1;}
		else {size_array[i] = divresult.quot;}
		temp--;
	}

///////////TESTOMRÅDE/////////////////////////////////////////



	//print2DArray(m,local_quot,b,rank); //print b på prosess 0
	//print1DArray(size, size_array, rank); //print size_array på prosess 0





// skriver ut max u for hver prosess
  	umax = 0.0;
  	for (j=0; j < local_quot; j++) {
    		for (i=0; i < m; i++) {
      			if (b[j][i] > umax) umax = b[j][i];
    		}
  	}
  	//printf (" umax = %e \n",umax);










//////////////////////////////////////////////////////////////

/*

  	transpose (bt,b,m);

  	for (i=0; i < m; i++) {
    		fstinv_(bt[i], &n, z, &nn);
  	}
  
  	for (j=0; j < m; j++) {
    		for (i=0; i < m; i++) {
      			bt[j][i] = bt[j][i]/(diag[i]+diag[j]);
    		}
  	}
  
  	for (i=0; i < m; i++) {
    		fst_(bt[i], &n, z, &nn);
  	}

  	transpose (b,bt,m);

  	for (j=0; j < m; j++) {
    		fstinv_(b[j], &n, z, &nn);
  	}

  	umax = 0.0;
  	for (j=0; j < m; j++) {
    		for (i=0; i < m; i++) {
      			if (b[j][i] > umax) umax = b[j][i];
    		}
  	}
  	printf (" umax = %e \n",umax);
	
	cout << "My rank is " << rank << endl;
*/  	
  	
  	MPI_Finalize();
  	return 0;
}
////////////////////////////////////////////////////////////
void print1DArray(int size, int array[], int rank){
	if (rank == 0) {
		for (int i = 0; i < size; i++){
			cout << array[i] << "\t";
		}
		cout << endl;
	}
}

void print2DArray(int m, int local_quot, Real **b, int rank){
	if (rank == 0) {
		cout << " Thread = " << rank << endl;
		for (int i=0;i<m;i++){
			for (int j=0; j < local_quot; j++) {
				cout << b[j][i] << "\t";
			}
			cout << "\n";
		}
	}
}

void transpose (Real **bt, Real **b, int m)
{
  int i, j;
  for (j=0; j < m; j++) {
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

