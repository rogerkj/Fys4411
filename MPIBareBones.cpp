/*
//
//          IsingMPI.cpp
//
//          Ising Thing
//
//          Roger Kjøde
//
*/


#include <mpi.h>
#include "lib.h"

#include <stdio.h>
#include <iostream>
#include <iomanip>

using namespace std;

//Used to store info about the parts of the calculation
struct parts {int element; bool isdone;};

//Used to send feedback to master proc
struct result {double sum;int element;int iam;};


//Entry point
int main(int argc, char* argv[])
{
 
//---  Changeable values  -----


	const int N = 10;

//------------------------------


	int numprocs, my_rank;
	MPI_Status status;

	//Init MPI
	MPI_Init (&argc, &argv);

	//Get number of procs and the rank of this proc
	MPI_Comm_size (MPI_COMM_WORLD, &numprocs);
	MPI_Comm_rank (MPI_COMM_WORLD, &my_rank);

	//Simple check for number of procs
	if (numprocs == 1) {
		printf("Number of nodes must be over 1!");
		MPI_Finalize ();
		return 0;
	}

//Master proc  Deligates tasks and collects feedback
	if(my_rank == 0) {

		//Used to store info of all the calculation parts
		parts* part = new parts[N];

		for(int t = 0;  t < N; t++) {

		  part[t].element = t;
		  part[t].isdone = false;
		  
	    
		}

		//If we have more processes than tasks
		int nextpart = numprocs - 1;
		if (N < numprocs) nextpart = N;

		//Send tasks to all procs
		for (int t2 = 1;  t2 <= nextpart ; t2++ ) {

			MPI_Send (&part[t2-1],sizeof(parts), MPI_BYTE,t2,100, MPI_COMM_WORLD);
		}

		//Used to store feedback
		result* answer = new result;

		//Used to tell procs to terminate.
		parts* donedummy = new parts;
		donedummy->isdone = true;

		int nrdone = 0;

		//Main feedback look
		do {
			//Wait for proc to send result
			MPI_Recv (answer,sizeof(result), MPI_BYTE, MPI_ANY_SOURCE, 100,MPI_COMM_WORLD, &status);

			//Store result
			
			part[answer->element].isdone = true;

		

			//Send proc new task or terminate proc if no more tasks available
			if(nextpart<N) {
				MPI_Send (&part[nextpart],sizeof(parts), MPI_BYTE,answer->iam, 100,MPI_COMM_WORLD);
			} else {
				MPI_Send (donedummy,sizeof(parts), MPI_BYTE,answer->iam, 100,MPI_COMM_WORLD);
			}

			//Index to next task incremented
			nextpart++;
			nrdone++;
			
		} while (nrdone < N);
		//Loop until all calculations are in

	
		printf("Done wone!");

		//All done!
		//All slave procs should now (or soon) be terminated.

	} else {
//Slave proc. Recives a task and replies an answer.

		//No use for me? Goodbye cruel world!
		if (N < my_rank) {
			MPI_Finalize();
			return 0;
		}

	

		//Used to store task and answer info
		parts* task = new parts;
		result* anser = new result;

		
		//Recive taks and send answer loop
		do {

			//Waits for new task.
			MPI_Recv (task, sizeof(parts), MPI_BYTE, 0 , 100,MPI_COMM_WORLD, &status);

			//Terminate proc. Not needed anymore
			if(task->isdone) break;

		
			anser->iam = my_rank;
			anser->element = task->element;

			//Send result to master proc 
			MPI_Send (anser, sizeof(result), MPI_BYTE, 0,100, MPI_COMM_WORLD);

		//Just loop
		} while(true);
	}

	//All done!
	MPI_Finalize ();

	return 0;
}
