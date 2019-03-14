// Thy Vu
// CSS 534 Autumn 2016
// Professor Munehiro Fukuda
// Programming Assignment 2

#include "mpi.h"
#include "omp.h"
#include <iostream>
#include "Timer.h"
#include <stdlib.h>   // atoi

int default_size = 100;  // the default system size
int defaultCellWidth = 8;
double c = 1.0;      // wave speed
double dt = 0.1;     // time quantum
double dd = 2.0;     // change in system

using namespace std;


int main( int argc, char *argv[] ) {
    
  int my_rank = 0;      // used by MPI
  int mpi_size = 1;     // used by MPI
  // verify arguments
  if ( argc < 4 ) {
    cerr << "usage: Wave2D size max_time interval" << endl;
    return -1;
  }
  int size = atoi( argv[1] );
  int max_time = atoi( argv[2] );
  int interval  = atoi( argv[3] );

  // get number of threads to run on each node
  int numThreads = 1;
  if ( argc == 5) {
    numThreads = atoi(argv[4]);
  }
  
  
  if ( size < 100 || max_time < 3 || interval < 0 ) {
    cerr << "usage: Wave2D size max_time interval" << endl;
    cerr << "       where size >= 100 && time >= 3 && interval >= 0" << endl;
    return -1;
  }

  MPI_Init( &argc, &argv ); // start MPI
  MPI_Comm_rank( MPI_COMM_WORLD, &my_rank );
  MPI_Comm_size( MPI_COMM_WORLD, &mpi_size );
  
  // create a simulation space
  double z[3][size][size];
  for ( int p = 0; p < 3; p++ ) 
    for ( int i = 0; i < size; i++ )
      for ( int j = 0; j < size; j++ )
	z[p][i][j] = 0.0; // no wave

  
  // defines each node's workload and work boundary
  int numColumnsPerNode = size / mpi_size;     
  int startColumn = my_rank * numColumnsPerNode;
  int endColumn = startColumn + numColumnsPerNode - 1;



  // start a timer
  Timer time;
  time.start( );

  // time = 0;
  // initialize the simulation space: calculate z[0][][]
  int weight = size / default_size;
  for( int i = 0; i < size; i++ ) {
    for( int j = 0; j < size; j++ ) {
      if( i > 40 * weight && i < 60 * weight  &&
	  j > 40 * weight && j < 60 * weight ) {
	z[0][i][j] = 20.0;
      } else {
	z[0][i][j] = 0.0;
      }
    }
  }

  // time = 1
  // calculate z[1][][] 
  // cells not on edge
  for (int i = 1; i < size-1; i++){
    for (int j = 1; j < size-1; j++){
        z[1][i][j] = z[0][i][j] + (c*c/2*dt*dt/dd/dd)*(z[0][i+1][j] + z[0][i-1][j]
                    + z[0][i][j+1] + z[0][i][j-1] - 4*z[0][i][j]);
    
    }
  }
  

  // simulate wave diffusion from time = 2 to max_time
  int now = 2;          // for z[] switching
  int nowMinus1 = 1;    // for z[] switching
  int nowMinus2 = 0;    // for z[] swtiching
  MPI_Status status;    // for MPI_Recv
  for ( int t = 2; t < max_time; t++ ) {      
  
    // Each node makes appropriate MPI calls to send its edge data 
    // and receive neighboring nodes' edge data
    if (my_rank == 0){  // master node
      // if there are multiple nodes, send my right edge 
      // to right neighbor and receive from the right neighbor.
      if (mpi_size > 1) {
        MPI_Send(z[nowMinus1][endColumn], size, MPI_DOUBLE, 1, 0, MPI_COMM_WORLD);
        MPI_Recv(z[nowMinus1][endColumn+1], size, MPI_DOUBLE,1, 0, MPI_COMM_WORLD, &status);
      }
      
    }else if (my_rank+1 == mpi_size){ // last node
      if (my_rank % 2 == 0){ // last even node
        // send my left, then receive from right neighbor
        MPI_Send(z[nowMinus1][startColumn], size, MPI_DOUBLE, my_rank-1, 0, MPI_COMM_WORLD);
        MPI_Recv(z[nowMinus1][startColumn-1], size, MPI_DOUBLE, my_rank-1, 0, MPI_COMM_WORLD, &status);
      
      } else { // last odd node
        // receive from left neighbor, then send my left edge
        MPI_Recv(z[nowMinus1][startColumn-1], size, MPI_DOUBLE, my_rank-1, 0, MPI_COMM_WORLD, &status);
        MPI_Send(z[nowMinus1][startColumn], size, MPI_DOUBLE, my_rank-1, 0, MPI_COMM_WORLD);
      }
  
    } else if (my_rank % 2 == 0){ // even-numbered node
      // send my left & my right
      MPI_Send(z[nowMinus1][startColumn], size, MPI_DOUBLE, my_rank-1, 0, MPI_COMM_WORLD);
      MPI_Send(z[nowMinus1][endColumn], size, MPI_DOUBLE, my_rank+1, 0, MPI_COMM_WORLD);
      
      // receive left & right from neighbors
      MPI_Recv(z[nowMinus1][startColumn-1], size, MPI_DOUBLE, my_rank-1, 0, MPI_COMM_WORLD, &status);
      MPI_Recv(z[nowMinus1][endColumn+1], size, MPI_DOUBLE, my_rank+1, 0, MPI_COMM_WORLD, &status);
   
    } else{   // odd-numbered node
      // receive left & right from neighbors
      MPI_Recv(z[nowMinus1][startColumn-1], size, MPI_DOUBLE, my_rank-1, 0, MPI_COMM_WORLD, &status);
      MPI_Recv(z[nowMinus1][endColumn+1], size, MPI_DOUBLE, my_rank+1, 0, MPI_COMM_WORLD, &status);
      
      // send my left & my right
      MPI_Send(z[nowMinus1][startColumn], size, MPI_DOUBLE, my_rank-1, 0, MPI_COMM_WORLD);
      MPI_Send(z[nowMinus1][endColumn], size, MPI_DOUBLE, my_rank+1, 0, MPI_COMM_WORLD);
    
    }
    
  
    // multithreaded schroedinger computation
    omp_set_num_threads(numThreads);
    #pragma omp parallel for 
    for (int i = startColumn; i <= endColumn ; i++){
        for (int j = 0; j < size-1; j++){
            if ( i != 0 && i != (size-1) && j != 0 && j != (size-1)) {
                z[now][i][j] = 2.0 * z[nowMinus1][i][j] - z[nowMinus2][i][j]+ (c*c*dt*dt/dd/dd) 
                            * (z[nowMinus1][i+1][j]+ z[nowMinus1][i-1][j] + z[nowMinus1][i][j+1] 
                            + z[nowMinus1][i][j-1]- 4.0*z[nowMinus1][i][j]);
            }
        }
    }
    
    // If it is time to print, slave nodes send computed data to master node,
    // master node receives and print.
    if (interval != 0 && (t % interval == 0 || t == max_time-1)) {
        if (my_rank != 0){
            // each slave node sends its computed results
            MPI_Send(z[now][startColumn], numColumnsPerNode * size, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
        } else {
            // master node receives
            for(int rank = 1; rank < mpi_size; rank++){
                MPI_Recv(z[now][rank * numColumnsPerNode], numColumnsPerNode * size, MPI_DOUBLE, rank, 0, MPI_COMM_WORLD, &status);
            }
           
            // master node prints
            cout<< t << endl;
            for (int i = 0; i < size; i++){
                for (int j = 0; j < size; j++){
                    cout<< z[now][i][j] << " ";
                }
                cout << endl;
            }
            cout<<endl;
            
        }
    }
    
    // switch z matrices
    int temp = now;
    now = nowMinus2;
    nowMinus2 = nowMinus1;
    nowMinus1 = temp;
        
   
  } // end of simulation

  // Print node rank and column range
  cerr << "rank["<< my_rank <<"]'s range = "<< startColumn<< " ~ "<< endColumn << endl;
  
  // finish the timer
  if (my_rank == 0)
      cerr << "Elapsed time = " << time.lap( ) << endl;
  
  MPI_Finalize();
  return 0;
}


