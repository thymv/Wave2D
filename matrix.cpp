#include "mpi.h"
#include <stdlib.h> // atoi
#include <iostream> // cerr
#include "Timer.h"

using namespace std;

void init( double *matrix, int size, char op ) {
  for ( int i = 0; i < size; i++ )
    for ( int j = 0; j < size; j++ )
      matrix[i * size + j]
	= ( op == '+' ) ? i + j :
	( ( op == '-' ) ? i - j : 0 );
}

void print( double *matrix, int size, char id ) {
  for ( int i = 0; i < size; i++ )
    for ( int j = 0; j < size; j++ )
      cout << id << "[" << i << "][" << j << "] = " << matrix[i * size + j] << endl;
}

void multiplication( double *a, double *b, double *c, int stripe, int size ) {
  for ( int k = 0; k < size; k++ )
    for ( int i = 0; i < stripe; i++ )
      for ( int j = 0; j < size; j++ )
	// c[i][k] += a[i][j] * b[j][k];
	c[i * size + k] += a[i * size + j] * b[j * size + k];
}

int main( int argc, char* argv[] ) {
  int my_rank = 0;            // used by MPI
  int mpi_size = 1;           // used by MPI
  int size = 400;             // array size
  bool print_option = false;  // print out c[] if it is true
  Timer timer;

  // variables verification
  if ( argc == 3 ) {
    if ( argv[2][0] == 'y' )
      print_option = true;
  }
  
  if ( argc == 2 || argc == 3 ) {
    size = atoi( argv[1] );
  }
  else {
    cerr << "usage:   matrix size [y|n]" << endl;
    cerr << "example: matrix 400   y" << endl;
    return -1;
  }

  MPI_Init( &argc, &argv ); // start MPI
  MPI_Comm_rank( MPI_COMM_WORLD, &my_rank );
  MPI_Comm_size( MPI_COMM_WORLD, &mpi_size );

  // matrix initialization
  double *a = new double[size * size];
  double *b = new double[size * size];
  double *c = new double[size * size];

  if ( my_rank == 0 ) { // master initializes all matrices
    init( a, size, '+' );
    init( b, size, '-' );
    init( c, size, '0' );

    // print initial values
    if ( false ) {
      print( a, size, 'a' );
      print( b, size, 'b' );
    }

    // start a timer
    timer.start( );
  }
  else {                // slavs zero-initializes all matrices
    init( a, size, '0' );
    init( b, size, '0' );
    init( c, size, '0' );
  }

  
  // broadcast the matrix size to all.
  MPI_Bcast( &size, 1, MPI_INT, 0, MPI_COMM_WORLD ); 

  int stripe = size / mpi_size;     // partitioned stripe. can assume divisible by mpi_size (eg 400/4)



  if (my_rank == 0) { //master sends [A] parts: 2nd stripe->rank1, 3rd->rk2, 4th->rnk3
  
    // master sends each partition of a[] to a different slave
    // master also sends b[] to all slaves
    
    for (int rank =1; rank < mpi_size; rank++){
        // MPI_send( void* message, int count, MPI_datatype, int dest, int tag, MPI_Comm comm)
        MPI_Send(a+rank*stripe*size , stripe * size, MPI_DOUBLE, rank, 0, MPI_COMM_WORLD);
        MPI_Send(b, size*size, MPI_DOUBLE, rank, 0, MPI_COMM_WORLD);
    }
  } else { 
    // MPI_Recv(void* message, int count, MPI_Datatype, int source, int tag, MPI_Comm, MPI_status* status)
    
    // slave receives partition of [a]
    // slave always put the received [A]'s stripe on its [A]'s top (easy: every slave computes from its A's top)
    MPI_Status status;
    MPI_Recv(a, stripe*size, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, &status);
    
    // slave receies [b]
    MPI_Recv(b, size*size, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, &status);
  
  }


  multiplication( a, b, c, stripe, size ); // all ranks should compute multiplication

  // master receives each partition of c[] from a different slave
  if (my_rank==0){
    for(int rank= 1; rank< mpi_size; rank++){ // receives from slave1->3
        MPI_Status status;
        MPI_Recv(c+rank*stripe*size, stripe *size, MPI_DOUBLE, rank , 0 , MPI_COMM_WORLD, &status);
    }
  
  }else{
    
    // slave sends its [c] computed portion
    MPI_Send(c, stripe*size, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD); // complete this tag=0?
  
  }

  if ( my_rank == 0 )
    // stop the timer
    cout << "elapsed time = " << timer.lap( ) << endl;

  // results
  if ( print_option && my_rank == 0 )
    print( c, size, 'c' );

  MPI_Finalize( ); // shut down MPI
}


  

  
