#include<iostream>
#include<cstdlib>
#include<cmath>
#include<pthread.h>
#include<fstream>
using namespace std;

long int num_threads = 2; //Number of threads
long int pattern_nodes,target_nodes,i,j,k;

//This holds the data that each thread requires:
struct thread_params{
	int part; //which part to multiply
	int M,N;
	long int **A;
	long int **B;
};

long int **mul;
long int **newP;
long int **H;
long int **targetMatrix;
long int **patternMatrix;
long int **Htranspose;
long int **X;

long int Test1();

float multiply_with_threads(long int** ,long int** ,int ,int);
void* matrix_multiply(void* params);

class AdjacencyMatrix
{
	private:
		long int n,p;
		bool *visited;

	public:
		long int matrixCount;
		long int **adj;	

		AdjacencyMatrix(long int n)
		{
			this->n = n;
			visited = new bool[n];
			adj = new long int* [n];
			for(long int i=0; i<n; i++)
			{
				adj[i] = new long int[n];
				for(long int j=0; j<n; j++)
				{
					adj[i][j]=0;
				}
			}
		}	

		//******************************* Function definition: void add_edge() *************************************//
		void add_edge(long int origin,long int destin)
		{
			if(origin>n || destin>n || origin<0 || origin<0)
			{
				cout<<"invalid edge"<<endl;
			}
			else
			{
				adj[origin][destin] = 1;
			}
		}//end of void add_edge

		//******************************* Function definition: void display() *************************************//
		void display(long int t,long int p,long int **dispMatrix)
		{
			for(long int i=0;i<p;i++)
			{
				for(long int j=0;j<t;j++)
					cout<<dispMatrix[i][j]<<"";
				cout<<endl;
			}

		}//end of void desplay

		//******************************* Function definition: void homoSearch() *************************************//		
		long int homoSearch(long int **H,long int **targetMatrix,long int **patternMatrix,long int target_nodes,long int pattern_nodes)
		{
			long int  rowIndex, colIndex; 	// rowIndex = i, columnIndex = j	
			long int x = pow(target_nodes, pattern_nodes); 
			long int breakCount=1;         //To come outof loop if match is not found!

			if(matrixCount == 1)
				rowIndex = -1;			//Build first leaf
			else rowIndex = pattern_nodes-1;	//H is already a leaf, build next leaf

			do
			{				
				if(breakCount == x)
				{
					break;	
				}

				if(rowIndex<(pattern_nodes-1))					
				{
					rowIndex = rowIndex + 1;

					for(colIndex=0; colIndex<target_nodes; colIndex++)
						if(colIndex == 0)
							H[rowIndex][colIndex] = 1;
						else
							H[rowIndex][colIndex] = 0;

				}
				else
				{
					while(rowIndex>=0 && H[rowIndex][target_nodes-1])
					{
						for(colIndex = 0; colIndex<target_nodes; colIndex++)
							H[rowIndex][colIndex] = 1; //true {restore i-th row}
					rowIndex = rowIndex - 1;	
					}//end of while	

					if (rowIndex == -1) //Matrix 'H' is null matrix				
					{
						for(long int i=0; i<pattern_nodes; i++)
						{
							for(long int j=0; j<target_nodes; j++)
								H[i][j]=0;						
						}
					}

					else{
						colIndex = 0;
						while(!H[rowIndex][colIndex])
						{
							colIndex = colIndex + 1;
						}
						H[rowIndex][colIndex] = 0; //false  //{unmap nodes i to node j}				
						H[rowIndex][colIndex+1] = 1; //true //{map node i to node j+!}
						breakCount++;
					}//end of else
				}//End of else

			}while(rowIndex!=(pattern_nodes-1) || !Test1() && rowIndex!=-1);

			if(breakCount == x)
			{
				cout<<"Pattern Not Found!!"<<endl;
			}
			else
			{
				cout<<"Pattern matrix found!!"<<endl;
				cout<<"Hence, The Homomorphism Matrix is:"<<endl;
				display(target_nodes,pattern_nodes,H);		
			}
		}//End of homoSearch function
};

//Test1 function:
long int Test1()
{
	AdjacencyMatrix adj(target_nodes);

	//Make mul matrix null before doing any multiplication.
	for(i=0; i<pattern_nodes; i++)
	{
		for(j=0; j<target_nodes; j++)
		{
			mul[i][j]=0;
			adj.matrixCount = 1;
		}
	}

	for(long int i=0; i<pattern_nodes; i++)
		for(long int j=0; j<target_nodes; j++)	
			Htranspose[j][i] = H[i][j];	
	
	//***************** Find H.T.Ht ***************//
	multiply_with_threads(H,targetMatrix,target_nodes,target_nodes);

	for(long int i=0; i<pattern_nodes; i++)
		for(long int j=0; j<target_nodes; j++)	
			X[i][j] = mul[i][j];

	for(i=0; i<pattern_nodes; i++)
	{
		for(j=0; j<target_nodes; j++)
		{
			mul[i][j]=0;
			adj.matrixCount = 1;
		}
	}

	multiply_with_threads(X,Htranspose,pattern_nodes,target_nodes);	

	//****************Check if newP is equal to old P**************//
	int flag = 0;	
	for(i=0; i<pattern_nodes; i++)
	{
		for(j=0; j<pattern_nodes; j++)
		{
			if(mul[i][j] != patternMatrix[i][j])
			{
				flag = 2;	
				break;
			}
			flag = 1;			
		}

		if(flag == 2)
			break;
	}

	if(flag == 1)
		return 1;
	else 
		return 0;
}//end of Test1


//function: multiply_with_threads
// Multiplies first and second using num_threads threads and puts the result in mat.
float multiply_with_threads(long int **A,long int **B, int P, int Q)
{
	int i,j,k;
	pthread_t* threads;  // pointer to a group of threads
	AdjacencyMatrix adj(target_nodes);

	threads = (pthread_t*) malloc(sizeof(pthread_t) * num_threads);

	// Create the threads
	for (i = 1; i < num_threads; i++)
	{

		//Reason for creating heap like below(line 305):
		//Local variables(here tparam is local to the thread) can go outof scope before thread starts executing.
		//hence create a heap and pass it to the thread!  
		//Do same thing for tparam->part =0;
		struct thread_params *tparam = (struct thread_params *) malloc(sizeof(struct thread_params));
		tparam->part = i;
		tparam->M = P;
		tparam->N = Q;
		tparam->A = A;
		tparam->B = B;
		
		if (pthread_create (&threads[i], NULL, matrix_multiply, tparam) != 0 )
		{
			cout<<"Error:pthread create error"<<endl;
			free(threads);
			exit(-1);
		}

	}

	struct thread_params *tparam = (struct thread_params *) malloc(sizeof(struct thread_params));
		tparam->part = 0;
		tparam->M = P;
		tparam->N = Q;
		tparam->A = A;
		tparam->B = B;

	matrix_multiply(tparam);

	// wait for all the threads to end
	for (i = 1; i < num_threads; i++)
		pthread_join (threads[i], NULL);
	free(threads);
}//multiply_with_threads

// A matrix is divided into 'part' equal partitions
// This is called by the thread with the 'part' number. 
void* matrix_multiply(void* params)
{
	AdjacencyMatrix adj(target_nodes);
	// Read the parameters
	struct thread_params *read_params;
	read_params = (thread_params*) params;
	int p = read_params->part;                                    // get the part number 

	int M = read_params->M;
	int N = read_params->N;
	long int **A = read_params->A;
	long int **B = read_params->B;

	int begin = (p * pattern_nodes)/num_threads;   // get the index for the beginning of this part
	int end = ((p+1) * pattern_nodes)/num_threads; // get the index for end of this part

	for(i=begin; i<end; i++)
		for(j=0; j<M; j++)
			for(k=0; k<N; k++)
				mul[i][j] = mul[i][j] || (A[i][k]*B[k][j]);				

			free(read_params);
	return 0;
}//end of matrix_multiply1

int main()
{
	long int i,j,t,max_edges, max_edges1, origin, destin,ex;
	char arr[20];
	ifstream myfile1;
	myfile1.open("example.txt");
	myfile1>>ex;
	cout<<"Enter number of examples:"<<ex<<endl;

for(t=0; t<ex; t++)
{	
	cout<<endl;
	myfile1>>arr;
	cout<<arr<<endl;;
	myfile1>>target_nodes;
	cout<<"Enter number of target nodes:"<<target_nodes<<endl;

	//**************************************Target Graph Creation***********************************//
	targetMatrix = new long int* [target_nodes];
	for(i=0; i<target_nodes; i++)
	{	
		targetMatrix[i] = new long int[target_nodes];
		for(j=0; j<target_nodes; j++)
		{
			targetMatrix[i][j]=0;
		}
	}

	AdjacencyMatrix am(target_nodes);
	max_edges = target_nodes*target_nodes;
	for(i=0;i<max_edges;i++)
	{
		myfile1>>origin>>destin;
		cout<<"Enter Edge:"<<origin<<" "<<destin<<endl;
		
		if((origin == -1) && (destin == -1))
			break;
		am.add_edge(origin-1,destin-1);
	}//end of for

	for(i=0;i<target_nodes;i++)
		for(j=0;j<target_nodes;j++)
			targetMatrix[i][j] = am.adj[i][j];

	cout<<"Target Matrix is:"<<endl;
	am.display(target_nodes,target_nodes,targetMatrix);	

	//***************************Pattern Graph Creation******************************//	
	myfile1>>pattern_nodes;
	cout<<"Enter number of pattern nodes:"<<pattern_nodes<<endl;
	
	patternMatrix = new long int* [pattern_nodes];
	for(i=0; i<pattern_nodes; i++)
	{
		patternMatrix[i] = new long int[pattern_nodes];
		for(j=0; j<pattern_nodes; j++)
		{
			patternMatrix[i][j]=0;
		}
	}
	
	//---mul, newP, Htranspose matrices initialization---//
	mul = new long int* [pattern_nodes];
	for(long int i=0; i<pattern_nodes; i++)
	mul[i] = new long int[target_nodes];

	X = new long int* [pattern_nodes];
	for(long int i=0; i<pattern_nodes; i++)
	X[i] = new long int[target_nodes];

	newP = new long int* [pattern_nodes];
	for(long i=0; i<pattern_nodes; i++)
	newP[i] = new long int[pattern_nodes];	

	Htranspose = new long int* [target_nodes];
	for(long int i=0; i<target_nodes; i++)
	Htranspose[i] = new long int[pattern_nodes];

	struct thread_params tp;

	tp.A = new long int* [pattern_nodes];
	for(long int i=0; i<pattern_nodes; i++)
	tp.A[i] = new long int[target_nodes];

	tp.B = new long int* [target_nodes];
	for(long int i=0; i<target_nodes; i++)
	tp.B[i] = new long int[target_nodes];

	AdjacencyMatrix am1(pattern_nodes);
	max_edges = pattern_nodes*(pattern_nodes-1);

	max_edges1 = pattern_nodes*pattern_nodes;
	for(i=0;i<max_edges1;i++)
	{
		myfile1>>origin>>destin;
		cout<<"Enter Edge:"<<origin<<" "<<destin<<endl;

		if((origin == -1) && (destin == -1))
			break;
		am1.add_edge(origin-1,destin-1);
	}//end of for

	for(i=0;i<pattern_nodes;i++)
		for(j=0;j<pattern_nodes;j++)
			patternMatrix[i][j] = am1.adj[i][j];
	
	cout<<"Pattern Matrix is:"<<endl;
	am1.display(pattern_nodes,pattern_nodes,patternMatrix);

	//*************************************Initial Homomorphism Matrix setup*************************/
	H = new long int* [pattern_nodes];
	for(i=0; i<pattern_nodes; i++)
	{
		H[i] = new long int[target_nodes];
		for(j=0; j<target_nodes; j++)
		{
			H[i][j]=1;
			am1.matrixCount = 1;
		}
	}

	am1.homoSearch(H, targetMatrix, patternMatrix, target_nodes, pattern_nodes);
	cout<<endl;
}//End of for

	myfile1.close();
	return 0;
}//End of main



