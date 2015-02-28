#include<mpi.h>
#include<stdio.h>
#include<string.h>
#include<time.h>

enum{
	MASTER_SEND = 0,
	SLAVE_SEND = 1,
};

//double estimate_g(double, double, long long int);
//void collect_results(double *); 

int main(int argc, char* argv[])
{
	MPI_Status status;
	
	
	float lower_bound = atof(argv[1]);
	float upper_bound = atof(argv[2]);
	long long int N = atof(argv[3]);
	
	//double sum;
	//sum = estimate_g(lower_bound, upper_bound, N);
	//printf("%lf is the final answer again\n",sum);
	
	int rank,cores;
	int source, dest;
	int average,extra;
	int i,j;
	int tag =0;
	double result=1,sum=0;
	double constant = (upper_bound - lower_bound)*(1/N)*(0.28209);
	
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &cores);
	
	average = N/(cores-1);
	extra = N%(cores-1);
	
	if(rank == 0)
	{
		for(dest=1; dest < cores; dest++)
		{
			if(extra!=0)
			{
				average++;
				MPI_Send(&average, 1, MPI_INT, dest, MASTER_SEND, MPI_COMM_WORLD);
				extra--;
				average--;
			}
			else
			{
			MPI_Send(&average, 1, MPI_INT, dest, MASTER_SEND, MPI_COMM_WORLD);
			}
		}
		printf("Sent all the data\n");
		for(source=1; source < cores; source++)
		{
			MPI_Recv(&result, 1, MPI_DOUBLE, source, SLAVE_SEND, MPI_COMM_WORLD, &status);
			sum = sum + result;
		}
		
		printf("Result is %lf from %d core\n",sum,rank);
			
	}
	
	else{
		for(j=1; j < cores; j++)
		{
		printf("Core %d has received the data\n",rank);
		MPI_Recv(&average, 1, MPI_INT, source, MASTER_SEND, MPI_COMM_WORLD, &status);
		int x = 3;
		for(i=0; i < average; i++)
		{
			result = result*3;
		}
		printf("%lf mid value\n",result);
		MPI_Send(&result, 1, MPI_DOUBLE, 0, SLAVE_SEND, MPI_COMM_WORLD);
		printf("Core %d has sent the data\n",rank);
		}
		
	}
    printf("If i see this message i am lucky\n");
    MPI_Finalize();                 // shuts down MPI
	//printf("this is stupid\n");
	return 0;
	
}	
/*
double estimate_g(double l_b, double u_b, long long int n)
{
	int rank,cores;
	int source, dest;
	int average,extra;
	int i,j;
	int tag =0;
	double result=0,sum=0;
	double constant = (u_b - l_b)*(1/n)*(0.28209);
	
	MPI_Status status;
	
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &cores);
	
	average = n/(cores-1);
	extra = n%(cores-1);
	
	if(rank == 0)
	{
		for(dest=1; dest < cores; dest++)
		{
			if(extra!=0)
			{
				average++;
				MPI_Send(&average, 1, MPI_INT, dest, MASTER_SEND, MPI_COMM_WORLD);
				extra--;
				average--;
			}
			else
			{
			MPI_Send(&average, 1, MPI_INT, dest, MASTER_SEND, MPI_COMM_WORLD);
			}
		}
		
		for(source=1; source < cores; source++)
		{
			MPI_Recv(&result, 1, MPI_DOUBLE, source, SLAVE_SEND, MPI_COMM_WORLD, &status);
			sum = sum + result;
		}
		
		printf("Result is %lf from %d core\n",sum,rank);
		return sum;
		
			
	}
	
	else{
	
		for(j=1; j < cores; j++)
		{
		MPI_Recv(&average, 1, MPI_INT, source, MASTER_SEND, MPI_COMM_WORLD, &status);
		//int x = 3;
		for(i=0; i < average; i++)
		{
			result++;
		}
		
		MPI_Send(&result, 1, MPI_DOUBLE, 0, SLAVE_SEND, MPI_COMM_WORLD);
		
		}
		
	}

}

*/


