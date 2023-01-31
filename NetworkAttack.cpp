#include <iostream>
#include <math.h>
#include <random>//For C++11 standard only.
#include <stdlib.h>
#include <time.h>//needed for randseed.
#include <fstream>
#include <stack>
#include <sstream>

using namespace std;

const int DEC_OUT = 12;
//decimal precision for console outputs.
 
const bool Debugging = false;

default_random_engine rng(7);
uniform_real_distribution<double> Pdist(0.0,1.0);

//*************************************************
//FORWARD DECLARATIONS.
//*************************************************
//---------NETWORK CONSTRUCTION----------------------
int CalculateCofKScaleFree(double NegGamma, int degmax, double CofK[]);
int CalculateCofKPoisson(double AvgK, int degmax,double CofK[]);
int PopulateDegreeList (int NodeN, int DegreeList[], double C_k[], int Kmax);
int DegreeListParity(int NodeN, int DegreeList[], int DegreeTotal, double C_k[], int Kmax);
int In_Out_Check (int NodeN, int DegreeList_i[], int DegreeTotal_i, int DegreeList_o[], int DegreeTotal_o, double C_k_i[],double C_k_o[], int KmaxIn, int KmaxOut);
int InitStartList(int Start[], int DegreeList[], int NodeN, int DegreeTotal);
int InitAdjacencyList(int Start[], int NodeN, int DegreeTotal, int Adj[]);
int ShuffleAList(int Adj[], int DegreeTotal);
int ShufflePerm2List(int Adj[], int DegreeTotal, int malpha);
int ShuffleUDAdjList(int Adj[], int DegreeTotal);
int BadEdgeRoutineWithPrune_Undirected(int Start[], int Adj[], int NodeN, int& DegreeTotal, int DegreeList[]);
int BadEdgeRoutineWithPrune_Directed(int Start_o[], int Start_i[], int Adj[], int NodeN, int& DegreeTotal, int DegreeList_o[], int DegreeList_i[]);

int BuildDirectedNetwork(int NodeN,int Kmax, char DistroType_i, char DistroType_o, int AddingModel, int run);
int BuildUndirectedNetwork(int NodeN,int Kmax,char DistroType, int AddingModel, double NetP_K[], int run);

//---------CENTRALITY MEASURE ALGORITHMS----------------
double eigenVector(int Start[], int NodeN, int Alist[],int DegreeTotal, double rightV[], bool leftVect);
int DynamicalImportance_Directed(int Start[], int Adjacency[], int DegreeTotal, int NodeN, double DynImp[]);
int DynamicalImportance_Undirected(int Start[], int Adjacency[], int DegreeTotal, int NodeN, double DynImp[]);
int BetCentContrib(int InStart[],int OutStart[], int NodeN, int Adjacency[], int DegreeTotal, double Bpart[], int root);
double *Betweenness(int InStart[],int OutStart[], int NodeN, int Adjacency[], int DegreeTotal, double Btotal[]);

//---------PERCOLATION FUNCTIONS-------------------
int FindClusterSize_U (int Start[], int Adjacency[], int del[], int visited[], int i, int ClusterSize);
int FindMaxCluster_U (int Start[], int Adjacency[], int del[], int NodeN);
int FindClusterSize_D (int Start[], int Adjacency[], int del[], int NodeN, int i);
//int FindMaxCluster_D (int Start[], int Adjacency[], int NodeN, int DegreeTotal);
void StrongConnect(int Start_o[], int Adjacency[], int visited[], int del[], int lowlink[], stack<int>& comp, int GSCC_vector[], int i, int& index, int add_index, int& gscc_length, int NodeN, int tempGSCC[], int inComp[]);
int GSCC(int Start_o[], int Adjacency[], int del[], int NodeN, int GSCC_vector[], int tempGSCC[]);

//--------UTILITIES-------------------------------
int EmpiricalPk(int DegreeList[], int DegreeTotal,double NetP_K[], int Kmax,int NodeN);
void printArray(int* array, int arraylen);
void printArray(double* array, int arraylen);
double PearsonCoefficient(int array1[], int array2[], int len);
double PearsonCoefficient(double array1[], double array2[], int len);

//---------NOISE ADDITION---------------------------------
void DeletingLinksUndirected (int StartTrue[], int StartFalse [], int AdjacencyTrue[], int AdjacencyFalse[], int Perm2[], int UnshuffledAdjacency[], double delta, int DegreeTotal, int NodeN);
void DeletingLinksDirected (int StartTrue[], int StartFalse [], int AdjacencyTrue[], int AdjacencyFalse[], int Perm2[], int DegreeList_i_False[], double delta, int DegreeTotal, int NodeN);
int PruneForFalse(int StartFalse[], int AdjacencyFalse[], int NodeN, int& DegreeTotalFalse, int DegreeListFalse[]);
void AddingEdgesModel1 (int Start[], int Adjacency[], int EdgesToAdd_1[], int EdgesToAdd_2[], double alpha, int NodeN, int EdgesToBeAdded, char BuildMode);
void AddingEdgesModel2_Undirected (int Start[], int Adjacency[], int UnshuffledAdjacency[], int EdgesToAdd_1[], int EdgesToAdd_2[], double alpha, int DegreeTotal, int EdgesToBeAdded);
void AddingEdgesModel2_Directed (int Start[], int Adjacency[], int UnshuffledAdjacency_i[], int UnshuffledAdjacency_o[], int EdgesToAdd_1[], int EdgesToAdd_2[], double alpha, int DegreeTotal, int EdgesToBeAdded);
int Partition_1(int EdgesToAdd_1[], int EdgesToAdd_2[], int p, int r);
void Quicksort_1(int EdgesToAdd_1[], int EdgesToAdd_2[], int p, int r);
int Partition_2(double EdgesToAdd_1[], int EdgesToAdd_2[], int p, int r);
void Quicksort_2(double EdgesToAdd_1[], int EdgesToAdd_2[], int p, int r);
void AddEdgesToStartAndAdjacency_Undirected (int FinalAdjacencyFalse[], int FinalStartFalse[], int AdjacencyFalse[], int StartFalse[], int DegreeListFalse[], int EdgesToBeAdded, int& DegreeTotalFalse, int NodeN, int EdgesToAdd_1[], int EdgesToAdd_2[]);
void AddEdgesToStartAndAdjacency_Directed (int FinalAdjacencyFalse[], int FinalStartFalse_o[], int FinalStartFalse_i[], int AdjacencyFalse[], int StartFalse_o[], int DegreeList_i_False[], int DegreeList_o_False[], int EdgesToBeAdded, int& DegreeTotalFalse, int NodeN, int EdgesToAdd_1[], int EdgesToAdd_2[]);

//-----------ATTACK ALGORITHMS--------------
int DegreeCentralityRanking(int DegreeList[], int NodeN, int DegCentrality[]);
void QuicksortCentralityRanking (double BTotal[], int NodeN, int BetweennessCentrality[]);
//---------RECALCULATIONS----------------
int DegreeListUpdate(int OrigDegreeList[], int Start[], int NodeN, int Adjacency[], int DegreeTotal, int del[], int NewDegList[]);
int BetCentContrib(int InStart[],int OutStart[], int NodeN, int Adjacency[], int DegreeTotal, double Bpart[], int root, int del[]);
double *Betweenness(int InStart[], int OutStart[], int NodeN, int Adjacency[], int DegreeTotal, double Btotal[], int del[]);
double eigenVector(int Start[], int NodeN, int Alist[],int DegreeTotal, double rightV[], bool leftVect,int del[]);
int DynamicalImportance_Directed(int Start[], int Adjacency[], int DegreeTotal, int NodeN, double DynImp[],int del[]);
int DynamicalImportance_Undirected(int Start[], int Adjacency[], int DegreeTotal, int NodeN, double DynImp[], int del[]);

//*************************************************
//----------------------------------------------------------


/*A program for the construction of directed and undirected complex networks, using a 'config model' given
 * degree distribution and node count, storing those networks in the form of an adjacency list.
 * 
 * Andrew Melchionna, Jesus Caloca, Shane Squires. June-Aug 2014 */
int main(int argc, char **argv){
	
	clock_t begin = clock();
	
	char BuildMode = 'u';
	char DistroType = 's';


	
	char DistroType_i = 's';
	char DistroType_o = 's';
	
	int AddingModel = 1;
	

	
	cout.precision(DEC_OUT);	
	srand(time(NULL));
	
	if (argc > 1){
		BuildMode = argv[1][1];
		//in "./Main.exe -u" picks out 'u' char.
	}
	if(Debugging == true)
		cout<<"BuildMode:"<< BuildMode <<endl;
	
	if(argc > 2){
		DistroType = argv[2][1];
	}
	if(Debugging == true)
		cout<<"Distribution Type(UD):"<< DistroType<< endl;	


	int Kmax = 1250;
	int NodeN = 2500;	
	
	int nruns = 50;

//for Undirected Networks
if (BuildMode == 'u'){
	
	double NetP_K[Kmax+1];
	double AvgScaleFreeArray[Kmax+1];
	
	
	for(int i=0;i<Kmax+1;i++)
	{AvgScaleFreeArray[i]=0.;}

	for (int run=1; run <= nruns;run++)
	{
		BuildUndirectedNetwork(NodeN, Kmax, DistroType, AddingModel, NetP_K, run);
  		for(int j=0;j<Kmax+1;j++)
		{AvgScaleFreeArray[j]+=NetP_K[j];}
		
		if(run % 500 == 0) {cout<<". ";
			if(run % 5000 == 0) cout<<endl;	}
		
	}
       FILE *f=fopen("dist.txt","w");
       for(int i=0;i<Kmax+1;i++)
       {
              AvgScaleFreeArray[i]/=((double)nruns);
              fprintf(f,"%15d%15.10f\n",i,AvgScaleFreeArray[i]);}
       fclose(f);	
}	

//Directed Networks
else if (BuildMode =='d'){	
	
	for (int run=1;run <= nruns;run++)
	{
	BuildDirectedNetwork(NodeN, Kmax, DistroType_i, DistroType_o, AddingModel, run);
		if(run % 500 == 0) {cout<<". ";
			if(run % 5000 == 0) cout<<endl;}
		
	}
}
	clock_t endd = clock();
	double difference = endd - begin;
	cout << "Elapsed Time: " << difference << endl;
	return 0;
}

//----------------------------------------------------------
//**********************************************************
//----------------------------------------------------------

//---------NETWORK CONSTRUCTION----------------------

/*TODO: Scale Free cumulative prob degre distribution*/
int CalculateCofKScaleFree(double NegGamma, int degmax, double CofK[]) {
	double ScaleFree[degmax+1];

	ScaleFree[0] = 0;
	CofK[0] =  ScaleFree[0];
	
	cout.precision(20);

	for (int k = 1;k <= degmax;k++)	{
		ScaleFree[k] = pow(k, NegGamma);
		}
	double ScaleFreeSum = 0.0;
	for (int i = 0; i <= degmax; i++) {
		ScaleFreeSum = ScaleFreeSum + ScaleFree[i];
	}
	double c = 1.0/ScaleFreeSum;
	for (int i = 0; i <= degmax; i++) {
		ScaleFree[i] = c * ScaleFree[i];
	}

	for (int k = 1;k <= degmax;k++)	{
	CofK[k] = CofK[k-1] + ScaleFree[k];
			if((1-CofK[k]) < 0 && k < degmax ){		
				cout<<"C(k) MAXED OUT AT "<< k<<endl; break; }

	}
	return 0;
}

/*C_KPOISSON Returns a cumulative Poisson distribution given
%an average degree count and a maximum degree count.
%INPUT:
%AvgK - The expectation/mean value that defines the Poisson distro.
%degmax - The maximum degree for which to calculate P(k), C(k) prob. distros.
%
%OUTPUT:
%CofK - The cumulative prob. distribution. Using matlab indices, CofK(1)
%   holds the cumulative prob, C(0).

%Pofx = ((AvgK)^k)*exp(-(AvgK))/factorial(k);
%amelchionna, jcaloca	*/
int CalculateCofKPoisson(double AvgK, int degmax,double CofK[]){
	double Poisson[degmax+1];	
	
	cout.precision(20);

//----------------
	Poisson[0] = 1;

	for (int k = 1;k <= degmax;k++)	{
		Poisson[k] = Poisson[k - 1] * AvgK/(double)k;
	}
	
	double C = 0.0;
	for (int i = 0; i <= degmax; i++) {
		C = C + Poisson[i];
		//cout << "Cequals " << C << endl;
	}
	
	Poisson[0] = Poisson[0]/C;
	CofK[0] = Poisson[0];
	for (int k = 1;k <= degmax;k++)	
	{
	
		Poisson[k] = Poisson[k]/C;
		CofK[k] = CofK[k-1] + Poisson[k];
		
		if((1-CofK[k]) < 0 && k < degmax ){		
			cout<<"C(k) MAXED OUT AT "<< k<<endl; break; }
		//FOR Both Debug and Non-debug modes, I would like to fail unambiguously when the rare combination of 
		//distribution and <double> values lead to an unreachable C(k)/ ultra high degree/improbable nodes.
			
	
	} 
	return 0;
}

/*Given any C(k) with a maximu Kmax, and a DegreeList of length NodeN, this function
 * distributes edge counts to the DegreeList using C(k)
 * */
int PopulateDegreeList (int NodeN, int DegreeList[], double C_k[], int Kmax){
	
	int DegreeTotal = 0;
	double r;
	
	for(int i = 0; i < NodeN;i++){		
	//1st node is at DegreeList[0]
	//last node at DegreeList[NodeN-1]

		r = Pdist(rng);
		
		for(int j=0;j<=Kmax;j++)	{			
			if(r<=C_k[j]){
				DegreeList[i] = j;
				DegreeTotal = DegreeTotal + j;
				break;
				}			
			}
			
	}//each node
	
	return DegreeTotal;
	}
	
int DegreeListParity(int NodeN, int DegreeList[], int DegreeTotal, double C_k[], int Kmax){
	
	if(Debugging==true)
		cout << "pre-parity:" << DegreeTotal<< "\n";
	double r;
	while (DegreeTotal % 2 != 0){
		int i = rand() % NodeN;//[0...NodeN-1] randomly
		DegreeTotal = DegreeTotal - DegreeList[i];
		
		r = Pdist(rng);
		
		for(int j=0;j<=Kmax;j++)	{			
			if(r<=C_k[j]){
				DegreeList[i] = j;
				DegreeTotal = DegreeTotal + j;
				break;
				}			
			}
	}
	if(Debugging == true)
		cout<< "Post-parity: " << DegreeTotal<< endl;
	
	return DegreeTotal;
	}

// CHECK FOR DEGREES IN = DEGREES OUT (DIRECTED ONLY)
int In_Out_Check (int NodeN, int DegreeList_i[], int DegreeTotal_i, int DegreeList_o[], int DegreeTotal_o, double C_k_i[],double C_k_o[], int KmaxIn, int KmaxOut)
	{

	int degree_sum_i_old = DegreeTotal_i;
	int degree_sum_o_old = DegreeTotal_o;
	int degree_sum_i_new;
	int degree_sum_o_new;

	while (degree_sum_i_old != degree_sum_o_old) {
		int r_node = rand() % NodeN;

		int d_i_new[1];
		int d_o_new[1];
		
		PopulateDegreeList(1, d_i_new, C_k_i, KmaxIn);		

		PopulateDegreeList(1, d_o_new, C_k_o, KmaxOut);	
		
		degree_sum_i_new = degree_sum_i_old + d_i_new[0] - DegreeList_i[r_node];
		degree_sum_o_new = degree_sum_o_old + d_o_new[0] - DegreeList_o[r_node];

		if (abs(degree_sum_i_new - degree_sum_o_new)
			< abs(degree_sum_i_old - degree_sum_o_old)) {
			degree_sum_i_old = degree_sum_i_new;
			degree_sum_o_old = degree_sum_o_new;
			DegreeList_i[r_node] = d_i_new[0];
			DegreeList_o[r_node] = d_o_new[0];
		}

		
		
	}
	
	return degree_sum_i_old;
	}

/*NOT USED. Given a probability r, and the C(k) cumulative prob. distribution, returns
 * a value of k corresponding to that r.*/
int BinProbSearch(double r, double* Cofk,int Kmax){
//TODO: Implement bin searching for the degree to assign. Will be more important for 
//SF networks.
return 0;	
}	

/*Given a list with the count of edges for each node, a 'start' list is constructed
with references the points in the adjacency list that have that node's edges*/
int InitStartList(int Start[], int DegreeList[], int NodeN, int DegreeTotal){
//	S = zeros(NodeN + 1,1,'uint16');

//S(1) = 1;
Start[0] = 0;

//for m = 2:NodeN
//Note: Start is N+1 long; the last element should index to the end of A list.
//The No. of edges for nth node is Start(n+1)-Start(n).
for(int s = 1;s<=NodeN;s++){
	Start[s] = Start[s-1] + DegreeList[s-1];
//    S(m) = S(m-1) + degree_list(m-1);	
}
//end
Start[NodeN] = DegreeTotal;
//S(NodeN + 1) = degree_sum + 1;
return 0;
}

/*Generate a List of 'stubs', which will correspond to the edges of the graph. 1st node 3 eges, 2nd node 1, we'll 
have: A[0-2] = 1, A[3] = 3 or
1	1	1	2	 */
int InitAdjacencyList(int Start[], int NodeN, int DegreeTotal, int Adj[]){
	//int i = 0;
	//for j = 2:NodeN
	int endind = 0;
	for(int i = 0;i < NodeN;i++){
		endind=i+1;
		for(int j = Start[i]; j < Start[endind];j++){
			Adj[j] = i;
		}
	}


return 0;
}

/*Shuffles the Perm List */
int ShuffleAList(int Adj[], int DegreeTotal){
	int SwapVal;
	int SwapIndex;	
	
	for(int i=0;i < DegreeTotal;i++){
		SwapIndex = (int) rand() % (DegreeTotal - i) + i;	

		SwapVal = Adj[SwapIndex];

		Adj[SwapIndex] = Adj[i];	

		Adj[i] = SwapVal;
	
	}

return 0;
}

/*Shuffles the Perm2 List for Noise Model (Deletion) */
int ShufflePerm2List(int Adj[], int EdgesTotal, int malpha){
	int SwapVal;
	int SwapIndex;

	for(int i=0;i < malpha ;i++){
		SwapIndex = (int) rand() % (EdgesTotal - i) + i;	

		SwapVal = Adj[SwapIndex];

		Adj[SwapIndex] = Adj[i];	

		Adj[i] = SwapVal;
	
	}

return 0;
}

/*Shuffles the initialized Adjacency List. */
int ShuffleUDAdjList(int Adj[], int DegreeTotal){
	
	if(DegreeTotal % 2 != 0){
	if (Debugging==true) cout<<"ShuffleUDAdjList called with odd numbered Adj list. No SHUFFLE."<<endl;
	return 1;
	}
	
	int Perm[DegreeTotal];
	int SwapVal = 0;
//obtain a randomized permutation of the edge indices.	
	for(int i=0;i<DegreeTotal;i++){
		Perm[i] = i;
	}	
	ShuffleAList(Perm, DegreeTotal);
	
//Do a pairwise swap so that nodes are connected in an undirected fashion.
	for(int i = 0;i < (DegreeTotal - 1);i=i+2){
		SwapVal = Adj[Perm[i+1]];
		Adj[Perm[i+1]]  = Adj[Perm[i]];
		Adj[Perm[i]] = SwapVal;
	}
	return 0;
}

/*This function will remove(deref & clear) self edges and double edges from the Adjacency list.
 * The number of bad edges will be returned.*/
int BadEdgeRoutineWithPrune_Undirected(int Start[], int Adj[], int NodeN, int& DegreeTotal, int DegreeList[]){
	//BadEdgeRoutineWithPrune(Start, Adjacency,NodeN,DegreeTotal);	
	
//	function [ BadCount, A ] = BadEdgeRoutineWithPrune(S, A)

int BadCount = 0;
int oldBadCount = 0;
bool jthIsBad = false;

//for i = 1:length(S)-1
for(int i = 0;i <= NodeN-1;i++){
	oldBadCount = BadCount;
	
//	for j = S(i):S(i+1)-1
	for(int j = Start[i];j < Start[i+1];j++){
//     Over each node's edge list.


		jthIsBad = false;
		
		if(i==Adj[j]){
		BadCount = BadCount+1;
		jthIsBad = true;
		if(Debugging == true)	
			cout<<"Self: "<< "("<<i<<", "<<j<<") => ("<<i<<", "<<Adj[j]<<") ."<<endl;
		} 
		else {
//        for k = j+1:S(i+1)-1
			/*Check this.*/
			for(int k = j+1; k<Start[i+1];k++){
				if(Adj[k]==Adj[j]){
               BadCount = BadCount+1;
               jthIsBad = true;      
      
               if(Debugging == true)
				cout<<"Double: "<< "("<<i<<", "<<k<<") = ("<<i<<", "<<Adj[k]<<") ."<<endl;

			   break;				
				}    
			}
		}//double edge check
		
		if(jthIsBad==false){
			Adj[j-BadCount] = Adj[j];	
		}//Dereference bad edges. 

	}//overnodes
Start[i] = Start[i] - oldBadCount;	
/*TODO: Does this neglect to alter the last element of Start?*/
}//over the start list.

Start[NodeN] = Start[NodeN] - BadCount;


	for (int j = DegreeTotal-(BadCount);j < DegreeTotal;j++){
		Adj[j] = -1;
	}

	for (int i = 0; i < NodeN; i++) {
		DegreeList[i] = (Start[i+1] - Start[i]);
	}

	DegreeTotal = DegreeTotal- BadCount;

	return BadCount;

}



int BadEdgeRoutineWithPrune_Directed(int Start_o[], int Start_i[], int Adj[], int NodeN, int& DegreeTotal, int DegreeList_o[], int DegreeList_i[]){
	//BadEdgeRoutineWithPrune(Start, Adjacency,NodeN,DegreeTotal);	
	
//	function [ BadCount, A ] = BadEdgeRoutineWithPrune(S, A)

int BadCount = 0;
int oldBadCount = 0;
bool jthIsBad = false;
int BadIn[NodeN+1];

for (int i = 0; i < (NodeN+1); i++) {
	BadIn[i] = 0;
}

//for i = 1:length(S)-1
for(int i = 0;i <= NodeN-1;i++){
	oldBadCount = BadCount;
	
//	for j = S(i):S(i+1)-1
	for(int j = Start_o[i];j < Start_o[i+1];j++){
//     Over each node's edge list.


		jthIsBad = false;
		
		if(i==Adj[j]){
		BadCount = BadCount+1;
		jthIsBad = true;
		BadIn[Adj[j]]++;
		if(Debugging == true)	
			cout<<"Self: "<< "("<<i<<", "<<j<<") => ("<<i<<", "<<Adj[j]<<") ."<<endl;
		} 
		else {
//        for k = j+1:S(i+1)-1
			/*Check this.*/
			for(int k = j+1; k<Start_o[i+1];k++){
				if(Adj[k]==Adj[j]){
               BadCount = BadCount+1;
               jthIsBad = true;  
			   BadIn[Adj[j]]++;
      
               if(Debugging == true)
				cout<<"Double: "<< "("<<i<<", "<<k<<") = ("<<i<<", "<<Adj[k]<<") ."<<endl;

			   break;				
				}    
			}
		}//double edge check
		
		if(jthIsBad==false){
			Adj[j-BadCount] = Adj[j];	
		}//Dereference bad edges. 

	}//overnodes
Start_o[i] = Start_o[i] - oldBadCount;	
/*TODO: Does this neglect to alter the last element of Start?*/
}//over the start list.

Start_o[NodeN] = Start_o[NodeN] - BadCount;


	for (int j = DegreeTotal-(BadCount);j < DegreeTotal;j++){
		Adj[j] = -1;
	}

	for (int i = 0; i < NodeN; i++) {
		DegreeList_o[i] = (Start_o[i+1] - Start_o[i]);
	}

	DegreeTotal = DegreeTotal- BadCount;

	int subtract_int = BadIn[0];
	for (int i = 1; i < NodeN; i++) {
		Start_i[i] -= subtract_int;
		subtract_int += BadIn[i];
	}

	Start_i[NodeN] -= subtract_int;

	for (int i = 0; i < NodeN; i++) {
		DegreeList_i[i] = (Start_i[i+1] - Start_i[i]);
	}

	return BadCount;

}

/**/
int BuildUndirectedNetwork(int NodeN,int Kmax,char DistroType, int AddingModel, double NetP_K[], int run){
	
	cout << "RUN " << run << endl;
	
	double C_k[Kmax+1];
	
	int DegreeList[NodeN];
	int DegreeListFalse[NodeN];
	int DegreeTotal;
	int DegreeTotalFalse;
	int BadCount;
	
	int Start[NodeN+1];	
	int StartFalse[NodeN + 1];
	
	double avgK = 4.0;
	
//	double Neggamma = -2.01421;//<k>=4 for kmax=500.
	double Neggamma = -2.06108;//<k>=4 for kmax=1250.
//	double Neggamma = -2.08569;//<k>=4 for kmax=2500.
//	double Neggamma = -2.10437;//<k>=4 for kmax=5000.


	double alpha = 0.0;
	double delta = 0.75;
	//Later this will have to be passed in.


//-----------------------



//-------------------------------


	if(DistroType == 's')
		CalculateCofKScaleFree(Neggamma,Kmax,C_k);
	if(DistroType == 'r')
		CalculateCofKPoisson(avgK, Kmax, C_k);
		
	if(Debugging==true){
		cout<<"C(k): " <<endl;
		printArray(C_k, Kmax+1);}
	
	if(Debugging==true)
		cout << "Populate DegreeList..."<<endl;
		
	DegreeTotal = PopulateDegreeList(NodeN, DegreeList,C_k,Kmax);
	//cout<<DegreeTotal<<endl;
	
	if(Debugging==true){
		cout<<"DegreeList: "<<endl;
		printArray(DegreeList, NodeN);	
	
		cout <<"Checking parity..."<<endl;	
		}
	
	DegreeTotal = DegreeListParity(NodeN, DegreeList, DegreeTotal, C_k,Kmax);	
	if(Debugging==true){
		printArray(DegreeList, NodeN);	
		
		cout << "DegreeTotal: " << endl;
		cout << DegreeTotal << endl;
		}	

	int Adjacency[DegreeTotal];
	int AdjacencyFalse[DegreeTotal];
	int Perm2[DegreeTotal/2];
	int UnshuffledAdjacency[DegreeTotal];
	

	InitStartList(Start, DegreeList, NodeN, DegreeTotal);
	
	
	if(Debugging==true){
		cout << "start:" <<endl;
		printArray(Start, NodeN + 1);
		}

		
	if(Debugging==true)	
		cout<<"Init AdjList..."<< endl;
	

	InitAdjacencyList(Start, NodeN, DegreeTotal, Adjacency);	

	
	if(Debugging==true) {
		printArray(Adjacency,DegreeTotal);
	
		cout<<"Shuffling AdjList: "<<endl; 
		}
	ShuffleUDAdjList(Adjacency, DegreeTotal);	
	
	if(Debugging==true){
		printArray(Adjacency, DegreeTotal);	
		cout<<"Removing Bad Edges..."<< endl;
	}



	BadCount = BadEdgeRoutineWithPrune_Undirected(Start, Adjacency,NodeN,DegreeTotal,DegreeList);
	if(Debugging==true){
		cout << "BadCount: " << endl;
		cout << BadCount << endl;
		cout << "Post-Prune DegreeTotal: " << endl;
		cout << DegreeTotal << endl;
		cout<<"Post-Prune StartList: "<< endl;
		printArray(Start, NodeN + 1);
		cout<<"Post-Prune AdjList: "<<endl;
		printArray(Adjacency, DegreeTotal);
	}
	cout << "BadCount: " << BadCount << endl;
	cout << "DegreeTotal: " << DegreeTotal << endl;
	
//---------CONCLUDES NET CONSTRUCTION---------


	
//Find the Constructed Network's Actual distribution.
	EmpiricalPk(DegreeList,DegreeTotal,NetP_K, Kmax, NodeN);
	
	if(Debugging==true){
		cout<<"Empirical Degree Distro,P(k)"<<endl;
		printArray(NetP_K,Kmax+1);		}
	
//Calculate the Average Degree.	
	double AvgDegreeEmpirical = 0;

	for (int i = 0; i < NodeN; i++) {
		AvgDegreeEmpirical = AvgDegreeEmpirical + DegreeList[i];
	}
	AvgDegreeEmpirical = AvgDegreeEmpirical / (double)NodeN;
	if(Debugging == true)cout << "AvgDegreeEmpirical: " << AvgDegreeEmpirical << endl;

	//NoiseModel
	
	for (int i = 0; i < NodeN; i++) {
		DegreeListFalse[i] = DegreeList[i];
	}

	DeletingLinksUndirected (Start, StartFalse, Adjacency, AdjacencyFalse, Perm2, UnshuffledAdjacency, delta, DegreeTotal, NodeN);
	DegreeTotalFalse = DegreeTotal;
	PruneForFalse(StartFalse, AdjacencyFalse, NodeN, DegreeTotalFalse, DegreeListFalse);


	//Adding Edges


	// We're actually adding degrees here
	double EdgesToBeRounded = DegreeTotal * alpha;
	int EdgesToBeAdded = (int) (EdgesToBeRounded + 0.5);
	int EdgesToAdd_1[EdgesToBeAdded];
	int EdgesToAdd_2[EdgesToBeAdded];
	int FinalAdjacencyFalse[DegreeTotalFalse + EdgesToBeAdded];
	int FinalStartFalse[NodeN+1];


	if (AddingModel == 1) {
		AddingEdgesModel1 (Start, Adjacency, EdgesToAdd_1, EdgesToAdd_2, alpha, NodeN, EdgesToBeAdded, 'u');
	}

	if (AddingModel == 2) {
		AddingEdgesModel2_Undirected (Start, Adjacency, UnshuffledAdjacency, EdgesToAdd_1, EdgesToAdd_2, alpha, DegreeTotal, EdgesToBeAdded);
	}

	


	Quicksort_1(EdgesToAdd_1, EdgesToAdd_2, 0, EdgesToBeAdded-1);
	


	AddEdgesToStartAndAdjacency_Undirected(FinalAdjacencyFalse, FinalStartFalse, AdjacencyFalse, StartFalse, DegreeListFalse, EdgesToBeAdded, DegreeTotalFalse, NodeN, EdgesToAdd_1, EdgesToAdd_2);


	
	
	
	if (Debugging == true) {
		cout << "AdjacencyFalse" << endl;
		printArray(AdjacencyFalse, DegreeTotalFalse - EdgesToBeAdded);
		cout << "StartFalse" << endl;
		printArray(StartFalse, NodeN + 1);

		cout << "EdgesToAdd_1" << endl;
		printArray(EdgesToAdd_1, EdgesToBeAdded);
		cout << "EdgesToAdd_2" << endl;
		printArray(EdgesToAdd_2, EdgesToBeAdded);

		cout << "FinalAdjacencyFalse" << endl;
		printArray(FinalAdjacencyFalse, DegreeTotalFalse);
		cout << "FinalStartFalse" << endl;
		printArray(FinalStartFalse, NodeN + 1);
	}
		
	//Centrality Measures

	//Percolation


	
	
	int DegCentrality[NodeN];
	for (int p = 0; p < NodeN; p++){
		DegCentrality[p] = p;
	}
	Quicksort_1(DegreeListFalse, DegCentrality, 0, NodeN-1);


	double Betw[NodeN];
	Betweenness(FinalStartFalse, FinalStartFalse , NodeN, FinalAdjacencyFalse, DegreeTotalFalse, Betw);

	int BetweennessCentrality[NodeN];
	QuicksortCentralityRanking (Betw, NodeN, BetweennessCentrality);	
	
	double DI[NodeN];
	DynamicalImportance_Undirected(FinalStartFalse,FinalAdjacencyFalse,DegreeTotalFalse,NodeN, DI);
	int DICentrality[NodeN];
	QuicksortCentralityRanking (DI, NodeN, DICentrality);

	int RandomDeletion[NodeN];
	for (int i = 0; i < NodeN; i++) {
		RandomDeletion[i] = i;
	}
	

	ShuffleAList(RandomDeletion, NodeN);
	
	


	
	//Percolation
	
	
	int del1[NodeN];
	int del2[NodeN];
	int del3[NodeN];
	int del4[NodeN];

	for (int i = 0; i < NodeN; i++) {
		del1[i] = 0;
		del2[i] = 0;
		del3[i] = 0;
		del4[i] = 0;
	}
	int MaxCluster1 = 1;
	int MaxCluster2 = 1;
	int MaxCluster3 = 1;
	int MaxCluster4 = 1;
	int kappa;
	
	double Common1;
	double Common2;
	double Common3;

	stringstream output;

	
	output << "50RUNS_" << "U_" << DistroType << "_" << 100*alpha << "_" << 100*delta << "_AM" << AddingModel << "_R" << run << ".txt";
	cout << output.str() << endl;

	ofstream a_file (output.str());
	
	MaxCluster1 = FindMaxCluster_U(Start, Adjacency, del1, NodeN);
	
	a_file <<  0.0 << "\t" << (double) MaxCluster1/NodeN << "\t";
	a_file << (double) MaxCluster1/NodeN << "\t" << (double) MaxCluster1/NodeN << "\t" << (double) MaxCluster1/NodeN << "\t";
	a_file << 0.0 << "\t"<< 0.0 << "\t"<< 0.0 << "\n";

	for (int i = NodeN - 1; i > 0; i--) {
		
		Common1 = 0.0;
		Common2 = 0.0;
		Common3 = 0.0;
		
		for (int k = NodeN - 1; k >= 0; k--) {
			if (del1[DegCentrality[k]] == 0) {
				del1[DegCentrality[k]] = 1;
				kappa = DegCentrality[k];
				break;
			}
		}
	
		MaxCluster1 = FindMaxCluster_U(Start, Adjacency, del1, NodeN);
		for (int already = 0; already < NodeN; already++) {
			if (del2[BetweennessCentrality[(NodeN-1) - already]] == 0) {
				del2[BetweennessCentrality[NodeN-1 - already]] = 1;
				break;

			}
		}
		
		MaxCluster2 = FindMaxCluster_U(Start, Adjacency, del2, NodeN);
		for (int already = 0; already < NodeN; already++) {
			if (del3[DICentrality[(NodeN-1) - already]] == 0) {
				del3[DICentrality[NodeN-1 - already]] = 1;
				break;

			}
		}
		MaxCluster3 = FindMaxCluster_U(Start, Adjacency, del3, NodeN);
		del4[RandomDeletion[i]] = 1;
		MaxCluster4 = FindMaxCluster_U(Start, Adjacency, del4, NodeN);
		
		for (int y = 0; y < NodeN; y++) {
			if (del2[y] == 1 && del3[y] == 1) {
				Common1 = Common1 + (double) 1/(NodeN - i);
			}
			if (del1[y] == 1 && del2[y] == 1) {
				Common2 = Common2 + (double) 1/(NodeN - i);
			}
			if (del1[y] == 1 && del3[y] == 1) {
				Common3 = Common3 + (double) 1/(NodeN - i);
			}
		}
				
		
		a_file << (double) (NodeN - i)/NodeN << "\t" << (double) MaxCluster1/NodeN << "\t";
		a_file << (double) MaxCluster2/NodeN << "\t" << (double) MaxCluster3/NodeN << "\t" << (double) MaxCluster4/NodeN << "\t";
		a_file << Common1 << "\t" << Common2 << "\t"<< Common3 << "\n";



		

		
		for (int j = FinalStartFalse[kappa]; j < FinalStartFalse[kappa+1]; j++) {
			DegreeListFalse[FinalAdjacencyFalse[j]] = DegreeListFalse[FinalAdjacencyFalse[j]] - 1;
		}
		
		
		Quicksort_1(DegreeListFalse, DegCentrality, 0, NodeN-1);
	
	


	Betweenness(FinalStartFalse, FinalStartFalse, NodeN, FinalAdjacencyFalse, DegreeTotalFalse, Betw, del2);
//	cout<<"Betweennesses, "<< i <<": "<<endl;
//	printArray(Betw,NodeN);

		
	QuicksortCentralityRanking (Betw, NodeN, BetweennessCentrality);
//	cout<<"Rankings, "<<i<<": "<<endl;
//	printArray(BetweennessCentrality, NodeN);
//

		
		
		DynamicalImportance_Directed(FinalStartFalse,FinalAdjacencyFalse,DegreeTotalFalse,NodeN, DI, del3);
		QuicksortCentralityRanking (DI, NodeN, DICentrality);
	 
	  	}

	a_file.close();
		


	/*int MaxCluster = FindMaxCluster_U(Start, Adjacency, del, NodeN);
	if(Debugging == true){
		cout << "MaxCluster is: ";
	}
	*/
	return 0;	
}

/**/
int BuildDirectedNetwork(int NodeN,int Kmax, char DistroType_i, char DistroType_o, int AddingModel, int run){
	
	double C_k_i[Kmax+1];
	double C_k_o[Kmax+1];
	
	int DegreeList_i[NodeN];
	int DegreeList_o[NodeN];
	int DegreeList_i_False[NodeN];
	int DegreeList_o_False[NodeN];
	int DegreeTotal_i;
	int DegreeTotal_o;
	int DegreeTotal;
	int DegreeTotalFalse;
	int BadCount;
	
	int Start_i[NodeN+1];
	int Start_o[NodeN+1];
	int Start_o_False[NodeN+1];
	
	double avgK = 4.0;

	//double Neggamma = -2.01471;

//	double Neggamma = -2.01421;//<k>=4 for kmax=500.
	double Neggamma = -2.06108;//<k>=4 for kmax=1250.
//	double Neggamma = -2.08569;//<k>=4 for kmax=2500.
//	double Neggamma = -2.10437;//<k>=4 for kmax=5000.


	double alpha = 1.0;
	double delta = 0.75;
	
//---------------------------------------

//Distribution Constructor
	if(DistroType_i == 'r')
		CalculateCofKPoisson(avgK,Kmax,C_k_i);
	if(DistroType_o == 'r')
		CalculateCofKPoisson(avgK,Kmax,C_k_o);
	if(DistroType_i =='s')
		CalculateCofKScaleFree(Neggamma,Kmax,C_k_i);
	if(DistroType_o == 's')
		CalculateCofKScaleFree(Neggamma,Kmax,C_k_o);
	
//-----------
	
	
	DegreeTotal_i = PopulateDegreeList(NodeN, DegreeList_i,C_k_i,Kmax);
	DegreeTotal_o = PopulateDegreeList(NodeN, DegreeList_o,C_k_o,Kmax);
	DegreeTotal = In_Out_Check (NodeN, DegreeList_i, DegreeTotal_i, DegreeList_o, DegreeTotal_o, C_k_i, C_k_o, Kmax, Kmax);
	int Adjacency[DegreeTotal];
	int AdjacencyFalse[DegreeTotal];
	int Perm2[DegreeTotal];
	int UnshuffledAdjacency_i[DegreeTotal];
	int UnshuffledAdjacency_o[DegreeTotal];
	
	if (Debugging == true) {
		cout << "DegreeList_i.."<<endl;	
		printArray(DegreeList_i,NodeN);
	
		cout << "DegreeList_o.."<<endl;	
		printArray(DegreeList_o, NodeN);
	
		cout << "Initial DegreeTotal: " << endl;
		cout << DegreeTotal << endl;
	}

	InitStartList(Start_o, DegreeList_o, NodeN, DegreeTotal);
	
	if (Debugging == true) {
		cout << "start_o:" << endl;	
		printArray(Start_o, NodeN + 1);
	}
	
	InitStartList(Start_i, DegreeList_i, NodeN, DegreeTotal);
	
	InitAdjacencyList(Start_i, NodeN, DegreeTotal, Adjacency);

	ShuffleAList(Adjacency, DegreeTotal);
	
	if (Debugging == true) {
		cout<<"Shuffled AdjList: "<<endl;		
		printArray(Adjacency, DegreeTotal);

		cout << "Start_i List" << endl;
		printArray(Start_i, NodeN + 1);
	}


	BadCount = BadEdgeRoutineWithPrune_Directed(Start_o, Start_i, Adjacency,NodeN,DegreeTotal,DegreeList_o, DegreeList_i);
	
	if (Debugging == true) {
	cout << "BadCount: " << endl;
	cout << BadCount << endl;
	cout << "Post-Prune DegreeTotal: " << endl;
	cout << DegreeTotal << endl;
	cout<<"Post-Prune Start_oList: "<<endl;
	printArray(Start_o, NodeN + 1);
	cout<<"Post-Prune Start_iList: "<<endl;
	printArray(Start_i, NodeN + 1);
	}
	
	if (Debugging == true) {
		cout<<"Post-Prune AdjList: "<<endl;	
		printArray(Adjacency, DegreeTotal);
	}
	
	cout << "Degree Total: " << DegreeTotal << endl;
	cout << "BadCount: " << BadCount << endl;
	
	//NoiseModel
	
	for (int i = 0; i < NodeN; i++) {
		DegreeList_i_False[i] = DegreeList_i[i];
		DegreeList_o_False[i] = DegreeList_o[i];
	}
	
	InitAdjacencyList(Start_i, NodeN, DegreeTotal, UnshuffledAdjacency_i);
	InitAdjacencyList(Start_o, NodeN, DegreeTotal, UnshuffledAdjacency_o);

	DeletingLinksDirected (Start_o, Start_o_False, Adjacency, AdjacencyFalse, Perm2, DegreeList_i_False, delta, DegreeTotal, NodeN);
	DegreeTotalFalse = DegreeTotal;
	PruneForFalse(Start_o_False, AdjacencyFalse, NodeN, DegreeTotalFalse, DegreeList_o_False);

	
	//Adding Edges

	double EdgesToBeRounded = DegreeTotal * alpha;
	int EdgesToBeAdded = (int) (EdgesToBeRounded + 0.5);
	int EdgesToAdd_1[EdgesToBeAdded];
	int EdgesToAdd_2[EdgesToBeAdded];
	int FinalAdjacencyFalse[DegreeTotalFalse + EdgesToBeAdded];
	int FinalStart_o_False[NodeN+1];
	int FinalStart_i_False[NodeN+1];

	if (AddingModel == 1) {
		AddingEdgesModel1 (Start_o, Adjacency, EdgesToAdd_1, EdgesToAdd_2, alpha, NodeN, EdgesToBeAdded, 'd');
	}
	if (AddingModel ==2) {
		AddingEdgesModel2_Directed (Start_o, Adjacency, UnshuffledAdjacency_i, UnshuffledAdjacency_o, EdgesToAdd_1, EdgesToAdd_2, alpha, DegreeTotal, EdgesToBeAdded);
	}


	Quicksort_1(EdgesToAdd_1, EdgesToAdd_2, 0, EdgesToBeAdded-1);


	AddEdgesToStartAndAdjacency_Directed(FinalAdjacencyFalse, FinalStart_o_False, FinalStart_i_False, AdjacencyFalse, Start_o_False, DegreeList_i_False, DegreeList_o_False, EdgesToBeAdded, DegreeTotalFalse, NodeN, EdgesToAdd_1, EdgesToAdd_2);

	if (Debugging == true) {
		
		cout << "Adjacency" << endl;
		printArray(Adjacency, DegreeTotal);
		cout << "Start" << endl;
		printArray(Start_o, NodeN + 1);
		
		cout << "AdjacencyFalse" << endl;
		printArray(AdjacencyFalse, DegreeTotal);
		cout << "StartFalse" << endl;
		printArray(Start_o_False, NodeN + 1);

		cout << "EdgesToAdd_1" << endl;
		printArray(EdgesToAdd_1, EdgesToBeAdded);
		cout << "EdgesToAdd_2" << endl;
		printArray(EdgesToAdd_2, EdgesToBeAdded);
		
		cout << "FinalAdjacencyFalse" << endl;
		printArray(FinalAdjacencyFalse, DegreeTotalFalse);
		cout << "FinalStartFalse" << endl;
		printArray(FinalStart_o_False, NodeN + 1);

		
	}
	

	if (Debugging == true) {
		cout << "DegreeListFalsei" << endl;
		printArray(DegreeList_i_False, NodeN); 
		cout<<"DegreeListFalseo"<<endl;
		printArray(DegreeList_o_False, NodeN);
	}


	int DegCentralityOut[NodeN];
	for (int p = 0; p < NodeN; p++){
		DegCentralityOut[p] = p;
	}
	Quicksort_1(DegreeList_o_False, DegCentralityOut, 0, NodeN-1);

	
	int DegCentralityIn[NodeN];
	for (int p = 0; p < NodeN; p++){
		DegCentralityIn[p] = p;
	}
	Quicksort_1(DegreeList_i_False, DegCentralityIn, 0, NodeN-1);



	double Betw[NodeN];
	Betweenness(FinalStart_i_False, FinalStart_o_False , NodeN, FinalAdjacencyFalse, DegreeTotalFalse, Betw);
	int BetweennessCentrality[NodeN];
	QuicksortCentralityRanking (Betw, NodeN, BetweennessCentrality);



	double DI[NodeN];
	DynamicalImportance_Directed(FinalStart_o_False,FinalAdjacencyFalse,DegreeTotalFalse,NodeN, DI);
	int DICentrality[NodeN];
	QuicksortCentralityRanking (DI, NodeN, DICentrality);
	
	int RandomDeletion[NodeN];
	for (int i = 0; i < NodeN; i++) {
		RandomDeletion[i] = i;
	}
	
	ShuffleAList(RandomDeletion, NodeN);
	
	


	//Percolation



	int del1[NodeN];
	int del2[NodeN];
	int del3[NodeN];
	int del4[NodeN];
	int del5[NodeN];


	for (int i = 0; i < NodeN; i++) {
		del1[i] = 0;
		del2[i] = 0;
		del3[i] = 0;
		del4[i] = 0;
		del5[i] = 0;
	}
	

	int GSCC1 = 1;
	int GSCC2 = 1;
	int GSCC3 = 1;
	int GSCC4 = 1;
	int GSCC5 = 1;
	

	int GSCC_vector1[NodeN];
	int GSCC_vector2[NodeN];
	int GSCC_vector3[NodeN];
	int GSCC_vector4[NodeN];
	int GSCC_vector5[NodeN];

	
	int tempGSCC1[NodeN];
	int tempGSCC2[NodeN];
	int tempGSCC3[NodeN];
	int tempGSCC4[NodeN];
	int tempGSCC5[NodeN];
	
	int OutCluster1;
	int OutCluster2;
	int OutCluster3;
	int OutCluster4;
	int OutCluster5;
	
	int region;
	int kappa1;
	int kappa2;
	
	double Common1;
	double Common2;
	double Common3;
	
	
	
	




	stringstream output;

	output << "DITest_D_" << DistroType_o << "_" << DistroType_i << "_" <<100*alpha << "_" << 100*delta << "_AM" << AddingModel << "_R" << run << ".txt";

	cout << output.str() << endl;

	ofstream a_file (output.str());
	
	GSCC1 = GSCC(Start_o, Adjacency, del1, NodeN, GSCC_vector1, tempGSCC1);
	OutCluster1 = FindClusterSize_D (Start_o, Adjacency, del1, NodeN, GSCC_vector1[0]);
	
	a_file << 0.0 << "\t" << (double) GSCC1/NodeN << "\t";
	a_file << (double) GSCC1/NodeN << "\t" << (double) GSCC1/NodeN << "\t" << (double) GSCC1/NodeN << "\t" << (double) GSCC1/NodeN << "\t";
	a_file <<(double) OutCluster1/NodeN << "\t" << (double) OutCluster1/NodeN << "\t" << (double) OutCluster1/NodeN;
	a_file << "\t" << (double) OutCluster1/NodeN << "\t" << (double) OutCluster1/NodeN << "\t" << 0.0 << "\t" << 0.0 << "\t" << 0.0 << "\n";
	



	for (int i = NodeN - 1; i > 0; i--) {
		Common1 = 0;
		Common2 = 0;
		Common3 = 0;
		
		for (int j = 0; j < NodeN; j++) {
			tempGSCC1[j] = NodeN;
			tempGSCC2[j] = NodeN;
			tempGSCC3[j] = NodeN;
			tempGSCC4[j] = NodeN;
			tempGSCC5[j] = NodeN;
		}

		for (int k = NodeN - 1; k >= 0; k--) {
			if (del1[DegCentralityOut[k]] == 0) {
				del1[DegCentralityOut[k]] = 1;
				kappa1 = DegCentralityOut[k];
				break;
			}
		}
		GSCC1 = GSCC(Start_o, Adjacency, del1, NodeN, GSCC_vector1, tempGSCC1);
		OutCluster1 = FindClusterSize_D (Start_o, Adjacency, del1, NodeN, GSCC_vector1[0]);
		for (int k = NodeN - 1; k >= 0; k--) {
			if (del2[DegCentralityIn[k]] == 0) {
				del2[DegCentralityIn[k]] = 1;
				kappa2 = DegCentralityIn[k];
				break;
			}
		}
		GSCC2 = GSCC(Start_o, Adjacency, del2, NodeN, GSCC_vector2, tempGSCC2);
		OutCluster2 = FindClusterSize_D (Start_o, Adjacency, del2, NodeN, GSCC_vector2[0]);
		for (int already = 0; already < NodeN; already++) {
			if (del3[BetweennessCentrality[(NodeN-1) - already]] == 0) {

				del3[BetweennessCentrality[NodeN-1 - already]] = 1;
				break;
			}
		}
		GSCC3 = GSCC(Start_o, Adjacency, del3, NodeN, GSCC_vector3, tempGSCC3);
		OutCluster3 = FindClusterSize_D (Start_o, Adjacency, del3, NodeN, GSCC_vector3[0]);
		for (int already = 0; already < NodeN; already++) {
			if (del4[DICentrality[(NodeN-1) - already]] == 0) {
				del4[DICentrality[NodeN-1 - already]] = 1;
				break;

			}
		}
		GSCC4 = GSCC(Start_o, Adjacency, del4, NodeN, GSCC_vector4, tempGSCC4);
		OutCluster4 = FindClusterSize_D (Start_o, Adjacency, del4, NodeN, GSCC_vector4[0]);
		del5[RandomDeletion[i]] = 1;
		GSCC5 = GSCC(Start_o, Adjacency, del5, NodeN, GSCC_vector5, tempGSCC5);
		OutCluster5 = FindClusterSize_D (Start_o, Adjacency, del5, NodeN, GSCC_vector5[0]);
		
		for (int y = 0; y < NodeN; y++) {
			if (del3[y] == 1 && del4[y] == 1) {
				Common1 = Common1 + (double) 1/(NodeN - i);
			}
			if (del1[y] == 1 && del3[y] == 1) {
				Common2 = Common2 + (double) 1/(NodeN - i);
			}
			if (del1[y] == 1 && del4[y] == 1) {
				Common3 = Common3 + (double) 1/(NodeN - i);
			}
		}

		a_file << (double) (NodeN - i)/NodeN << "\t" << (double) GSCC1/NodeN << "\t";
		a_file << (double) GSCC2/NodeN << "\t" << (double) GSCC3/NodeN << "\t" << (double) GSCC4/NodeN << "\t" << (double) GSCC5/NodeN << "\t";
		a_file <<(double) OutCluster1/NodeN << "\t" << (double) OutCluster2/NodeN << "\t" << (double) OutCluster3/NodeN;
		a_file << "\t" << (double) OutCluster4/NodeN << "\t" << (double) OutCluster5/NodeN << "\t";
		a_file << Common1 << "\t"<< Common2 << "\t" << Common3 << "\n";
		

		
		for (int j = 0; j < DegreeTotalFalse; j++) {
			if (FinalAdjacencyFalse[j] == kappa1) {
				for (int k = 0; k < NodeN; k++) {
					if (FinalStart_o_False[k] > j) {
						region = k - 1;
						break;
					}
				}
			DegreeList_o_False[region]--;
			}
		}
		Quicksort_1(DegreeList_o_False, DegCentralityOut, 0, NodeN-1);
		
		
				
		
		for (int j = FinalStart_o_False[kappa2]; j < FinalStart_o_False[(kappa2)+1]; j++) {
			DegreeList_i_False[FinalAdjacencyFalse[j]] = DegreeList_i_False[FinalAdjacencyFalse[j]] - 1;
		}
		
		Quicksort_1(DegreeList_i_False, DegCentralityIn, 0, NodeN-1);
		

		Betweenness(FinalStart_i_False, FinalStart_o_False , NodeN, FinalAdjacencyFalse, DegreeTotalFalse, Betw, del3);
		QuicksortCentralityRanking (Betw, NodeN, BetweennessCentrality);
		
		
		
		
		DynamicalImportance_Directed(FinalStart_o_False,FinalAdjacencyFalse,DegreeTotalFalse,NodeN, DI, del4);
		QuicksortCentralityRanking (DI, NodeN, DICentrality);


	
	}
	a_file.close();
	
	return 0;
}

//******************************************************
//---------CENTRALITY MEASURE ALGORITHMS----------------

/*From the Start and Adjacency lists, determine the eigenvector corresponding to the largest
 * eigenvalue of this matrix. The 'eigVector' passed will have the vector filled upon return. 
The largest eigenvalue is also returned as a double by the function. 
 * FOR DIRECTED NETWORKS: Pass-in the out_Start list.(From the derivation of the Dynamical 
 * Importance approx., Restrepo,et al 2006: "We consider a network as a directed graph with N nodes,
and we associate to it a N x N adjacency matrix whose
elements Aij are positive if there is a link going from node i
to node j with i != j and zero otherwise (Aii = 0)."		 */
double eigenVector(int Start[], int NodeN, int Alist[], int DegreeTotal,double eigVector[],bool leftVect){
	
	//Prepare a pair of vectors of length N, for iteration. Let the iteration begin with 
	//a random vector, which is then normalized.
	double newV[NodeN];	
	double oldV[NodeN];
	double NormSum = 0.0;
	
	int a = 0;
	int power = 0;	
	
	double nsL=0;//New and Old Shifted Lambdas(\lambda + 1 = nsL)
	double osL=0;
	double Lambda = 0;
	
	
	//init a blank V_n+1, let V_0 be a random vector, which is normalized.
	for(int i = 0;i<NodeN;i++) {
		newV[i] = 0.0;	
		oldV[i] = Pdist(rng);		
		NormSum += oldV[i]*oldV[i];
	}
		
	NormSum = sqrt(NormSum);		
	for(int i = 0;i < NodeN;i++) oldV[i] /= NormSum;	
	osL = NormSum;
	
	if(Debugging==true){
	cout<<"V_"<< 0 << endl;
	printArray(oldV,NodeN);	}

do{//while lambda has not converged to within 1e-5.
	
	power++;
	//Multiply the Adjacency matrix to oldV to give newV.
	for(int i =0; i < NodeN;i++){//V_n+1 = V_n + A*V_n
		newV[i] = oldV[i];
	}
	for(int i=0; i < NodeN; i++){
		for(int j=Start[i];j<Start[i+1];j++){
			a = Alist[j];
				
			if(leftVect==false)newV[i] += oldV[a];
				else newV[a] += oldV[i];
		//Per notes: The left eigenvector of A is equal to the right 
		//	eigenvector of A^T(transpose). The switch of indices uses the transpose.  
		}	
	}	
	
	//Calculate the Norm(root(Sum(V)) for n+1.
	for(int i =0; i<NodeN;i++) nsL += newV[i]*newV[i];
	nsL = sqrt(nsL);
	
	//normalize the new V_n+1
	for(int i =0; i < NodeN;i++) oldV[i] = newV[i] / nsL;
	//cout << power << "   " << nsL << "   " << newV[0] << "   " << nsL-osL << endl;
	
	//cout<<"nsL = " << nsL <<". V_"<< power <<": "<< endl;
	//printArray(oldV,NodeN);	
	
	//Check for convergence of the eigenvalue; The magnitude of the V_n+1 we obtain.	
	if(fabs(nsL - osL) < 1e-5) break;
	else {	osL = nsL;
			nsL = 0;		}
	
}while(true);

	Lambda = nsL - 1.0;	
	if(Debugging==true)	cout<<"Found EigenValue- "<<Lambda<<endl;
	
	if(Debugging==true){ cout<<"Found EigenVector: \n";printArray(oldV,NodeN);}
	for(int p =0;p<NodeN;p++) eigVector[p] = oldV[p];	
	
	return Lambda;
	
}

/*Using the approximation given in Restrepo, et al. 2006:
 * I(k) = u_k*v_k/(u^T * v ) The DI of the kth node is approximated by taking the product
 * of the kth elements of the left and right eigenVectors and normalizing them by the dot product
 * of those vectors. The given array, DynImp[] will be filled with these values.  */
int DynamicalImportance_Directed(int Start[], int Adjacency[], int DegreeTotal, int NodeN, double DynImp[]){
	//Take the dot product of the left/right eigenvectors. This is the normalization factor.
	
	double lVect[NodeN];
	double rVect[NodeN];
	eigenVector(Start, NodeN, Adjacency, DegreeTotal, lVect, true);
	eigenVector(Start, NodeN, Adjacency, DegreeTotal, rVect, false);
	
	double NormFactor= 0.0;	
	for(int i = 0;i<NodeN;i++){	NormFactor += lVect[i] * rVect[i];	}
	//Note: Should be 1.0 for the undirected case; as the eigVector should be normalized before call.
	
	//Into DynImp, calculate the approximate DI for every node, as given above.
	for(int i = 0;i<NodeN;i++){ DynImp[i] = lVect[i] * rVect[i] / NormFactor;	}	
	
	return 0;
}


int DynamicalImportance_Undirected(int Start[], int Adjacency[], int DegreeTotal, int NodeN, double DynImp[]){
	
	double eVect[NodeN];	
	eigenVector(Start, NodeN, Adjacency, DegreeTotal, eVect, false);
	
	double NormFactor= 0.0;	
	for(int i = 0;i<NodeN;i++){	
		NormFactor += eVect[i] * eVect[i];	
		}
	
	//Into DynImp, calculate the approximate DI for every node, as given above.
	for(int i = 0;i<NodeN;i++){
	DynImp[i] = eVect[i] * eVect[i] / NormFactor;
	}	

	return 0;
}

/*For the given root node, will calculate all of the shortest paths,
 * and use this list to calculate node's betweenness.*/
int BetCentContrib(int InStart[],int OutStart[], int NodeN, int Adjacency[], int DegreeTotal, double Bpart[], int root){

// q: a queue of nodes to visit
// qs: start of used portion of q
// qe: end of used portion of q
	int Q[NodeN];
	int qs=0;int qe=0;
	
//distance to a node
	int DtoRoot[NodeN];
	
//number of paths of minimum distance to each node
	int PathsFromRoot[NodeN];
  
//!Initialize queue to contain only root node, distances to be max 
//!possible, and betweennesses to be 0.
 for(int n = 0;n < NodeN;n++){
	 Q[n]=NodeN+1;
	 DtoRoot[n] = NodeN+1;
	 Bpart[n] = 0;
	 PathsFromRoot[n]=0;
 }
 Q[0]=root;
 DtoRoot[root] = 0;
 PathsFromRoot[root] = 1; 
 
int parents[DegreeTotal];
for(int p = 0; p< DegreeTotal;p++){
	parents[p] = NodeN+1;}

//A Breadth First Search for the shortest paths to all other nodes that begin with 
//the root given. 
 int i, j, IndIncr;
 do 
 {
	i = Q[qs];
	qs = qs + 1;
		//Directed Version; branch out uses InStart, parent check uses OutStart
		for(int e = OutStart[i];e<OutStart[i+1];e++){	
			j = Adjacency[e];
			if(DtoRoot[j] >= (DtoRoot[i]+1)){//condition for a shortest path, possibly non-unique
				if(DtoRoot[j] > DtoRoot[i]+1){//condition for unexplored node
					qe=qe+1;
					Q[qe] = j;
					DtoRoot[j] = DtoRoot[i] + 1;
				}
				PathsFromRoot[j] = PathsFromRoot[j] + PathsFromRoot[i];
				
			//For that node's particular parent list, only add if no shortest path
			//parent occupies that position already.
					IndIncr = 0;
					while(parents[InStart[j]+IndIncr] < NodeN+1 ){IndIncr++; }
					parents[InStart[j]+IndIncr] = i;
				}			
		}
//cout<<"Paths from root:"<<endl;
//printArray(PathsFromRoot, NodeN);		
//cout<<"Distances to r("<<root<<"):"<<endl;
//printArray(DtoRoot, NodeN);	
//cout<<"Queue("<<qs<<"..."<<qe<<"):"<<endl;
//printArray(Q, NodeN);
//cout<<"Parents:"<<endl;
//printArray(parents,DegreeTotal);
//cout<<"-->Next:"<<Q[qs]<<endl;
	 
}while (qs <= qe);

//------------------- 

//Betweenness, given num paths to root, distances to root, and paths via parent trace-back;
//calculates Betweenness Centrality Contribution for all nodes != to the root.
int k;
for(int Queuei = qe; Queuei>= 0; Queuei--){
	IndIncr = 0;

	k=Q[Queuei];

		while(parents[InStart[k]+IndIncr]<NodeN+1 && InStart[k]+IndIncr < InStart[k+1]){
	//pull parents from the part of parents list that is filled, but not going over into the
	//next node's parents.
		j = parents[InStart[k]+IndIncr];
		Bpart[j] = Bpart[j]+ (1.0 + Bpart[k] )*(double)PathsFromRoot[j] / (double)PathsFromRoot[k];
		IndIncr++;
	}		
}  
 Bpart[root] = 0; 
 
 if(Debugging==true)
	printArray(Bpart,NodeN);
 //working backward from the end of the queue. Algorithm per S.S.	
//INITIALIZE: bb(i)=0 for all i
//for each i (going backwards in the queue)
//for each parent of i,
//       bb(parent)=bb(parent)+(1.+bb(i))*((float)np(parent))/((float)np(i))
//AFTERWARDS: bb(root node)=0.
 return 0;

}

/*For both directed and undirected networks, calculates the betweenness of each node by counting the number of times
 * each nodes occurs in the shortest paths between all pairs(weighting when there is more than 1 s.p.). 
 * IN: In and Out Start lists for Directed. FOR UNDIRECTED, SEND 'START' for both inputs. 
 * 'Btotal' will have Betweenness values upon return.  * */
 double *Betweenness(int InStart[], int OutStart[], int NodeN, int Adjacency[], int DegreeTotal, double Btotal[]){
	double Bpart[NodeN];
	//every Node will run BetCentContrib, using shortest paths from that node to all
	//other nodes to determine the contribution of paths from that node to the betweenness
	//of every other node.	

	for(int i=0;i<NodeN;i++) Btotal[i]=0.;
	for(int root = 0;root < NodeN;root++)
	{
		for(int i = 0;i< NodeN;i++){
		Bpart[i] = 0.0;
		}

		if(Debugging==true)	cout<<"Contrib, Node "<<root<<endl;
		BetCentContrib(InStart, OutStart, NodeN, Adjacency, DegreeTotal, Bpart, root);

		for(int j = 0;j<NodeN;j++){
		Btotal[j] = Btotal[j] + Bpart[j];
		}
		
		if(root % 1000 == 0){
				cout<<"* ";
					if(root % 15000==0){cout<<endl;}
						}
	}
	 
	if(Debugging==true)
		{cout<<"Betweenness:"<<endl; printArray(Btotal, NodeN);}
	
	return Btotal;
}
//*****************************************************
//---------PERCOLATION FUNCTIONS-------------------

//MaxClusterSize for undirected
int FindClusterSize_U (int Start[], int Adjacency[], int del[], int visited[], int i, int ClusterSize) {
	int k;
	visited[i] = 1;
	ClusterSize = ClusterSize + 1;

	for (int j = Start[i]; j < Start[i + 1]; j++) {
		if (del[Adjacency[j]] == 0 && visited[Adjacency[j]] == 0) {
			k = Adjacency[j];
			ClusterSize = FindClusterSize_U (Start, Adjacency, del, visited, k, ClusterSize);
		}
	}
	return ClusterSize;
}

int FindMaxCluster_U (int Start[], int Adjacency[], int del[], int NodeN){

	int visited[NodeN];
	
	//Initialize del list
	for (int i = 0; i < NodeN; i++) {
		visited[i] = 0;		}
	
	int ClusterSize = 0;
	int MaxCluster = 0;
	for (int i = 0; i < NodeN; i++) {
		if (visited[i] == 0 && del[i] == 0) {
			ClusterSize = FindClusterSize_U (Start, Adjacency, del, visited, i, 0);
		}
		if (ClusterSize > MaxCluster) {
			MaxCluster = ClusterSize;
		}
		if(Debugging==true){
			cout << "Visited:" << endl;
			for (i = 0; i < NodeN; i++) {
			cout << visited[i] << endl;
			}
		}
		
	}
	return MaxCluster;
}

//MaxClusterSize For Directed
int FindClusterSize_D (int Start[], int Adjacency[], int del[], int NodeN, int i) {
	int visited[NodeN];
	int queue[NodeN + 1];
	int t;
	
		//Initialize del list
	
	for (int i = 0; i < NodeN; i++) {
		visited[i] = 0;
	}


	for (int k = 0; k < NodeN + 1; k++) {
		queue[k] = 0;
	}
	
	queue[0] = (i+1);
	visited[i] = 1;
	int q_index = 1;
	int ClusterSize = 1;
	while (queue[0] != 0) {
		t = queue[0] - 1;
		q_index = q_index - 1;
		for (int z = 0; z < NodeN; z++) {
			queue[z] = queue[z+1];
		}
		for (int j = Start[t]; j < Start[t + 1]; j++) {
			if (del[Adjacency[j]] == 0 && visited[Adjacency[j]] == 0) {
				visited[Adjacency[j]] = 1;
				ClusterSize++;
				queue[q_index] = (Adjacency[j] + 1);
				q_index++;
			}
		}
	}
	return ClusterSize;
}

/*int FindMaxCluster_D (int Start[], int Adjacency[], int NodeN, int DegreeTotal) {

	int visited[NodeN];
	int del [NodeN];
	int queue[NodeN + 1];
	
		//Initialize del list
	for (int i = 0; i < NodeN; i++) {
		del[i] = 0;
	}
	
	for (int i = 0; i < NodeN; i++) {
		visited[i] = 0;
	}
	int ClusterSize = 0;
	int MaxCluster = 0;

	for (int i = 0; i < NodeN; i++) {
		if (visited[i] == 0 && del[i] == 0) {
			ClusterSize = FindClusterSize_D (Start, Adjacency, del, visited, queue, NodeN, i, it_index);
			if (Debugging==true) cout << "Intermediate MaxCluster: " << ClusterSize << endl;

			if (ClusterSize > MaxCluster) {
				MaxCluster = ClusterSize;
			}
		}
	}

	return MaxCluster;
}*/


//Giant Strongly Connected Component
void StrongConnect(int Start_o[], int Adjacency[], int visited[], int del[], int lowlink[], stack<int>& comp, int GSCC_vector[], int i, int& index, int add_index, int& gscc_length, int NodeN, int tempGSCC[], int inComp[]) {
	//stack<int> temp;
	visited[i] = index;
	lowlink[i] = index;
	index++;
	comp.push(i);
	inComp[i] = 1;
	
	
	for (int j = Start_o[i]; j < Start_o[i + 1]; j++) {
/*		bool EdgeInStack;

		while (!comp.empty() && comp.top() != Adjacency[j]) {
			temp.push(comp.top());
			comp.pop();
		}

		while (!temp.empty()) {
			comp.push(temp.top());
			temp.pop();
		}

		if (comp.empty()) {
				EdgeInStack = false;
			}
		if (comp.top() == Adjacency[j]) {
			EdgeInStack = true;
		}*/
		
		if (visited[Adjacency[j]]==0 && del[Adjacency[j]] == 0) {
			StrongConnect(Start_o, Adjacency, visited, del, lowlink, comp, GSCC_vector, Adjacency[j], index, add_index, gscc_length, NodeN, tempGSCC, inComp);
			if (lowlink[i] > lowlink[Adjacency[j]] && inComp[Adjacency[j]] == 1) {
				lowlink[i] = lowlink[Adjacency[j]];
			}
		}
		else if (del[Adjacency[j]]==0) {
			if (lowlink[i] > visited[Adjacency[j]] && inComp[Adjacency[j]] == 1) {
				lowlink[i] = visited[Adjacency[j]];
			}
		}
	}

	if (lowlink[i] == visited[i]) {
		

		while (comp.top() != i) {
			tempGSCC[add_index] = comp.top();
			inComp[comp.top()] = 0;
			comp.pop();
			add_index++;
				
		}
		
		tempGSCC[add_index] = i;
		inComp[i] = 0;
		comp.pop();
		if(Debugging==true)cout << "here's tempgscc" << endl;
		for (int k = 0; k <= add_index; k++) {
			if(Debugging==true)cout<<tempGSCC[k]<<endl;
		}

		if (add_index  >= gscc_length) {
			for (int k = 0; k <= add_index; k++) {
				GSCC_vector[k] = tempGSCC[k];
			}
			gscc_length = add_index;
		}
	}
}


int GSCC(int Start_o[], int Adjacency[], int del[], int NodeN, int GSCC_vector[], int tempGSCC[]) {

	int index = 1;
	int gscc_length = 0;
	int visited[NodeN];
	int lowlink[NodeN];
	stack<int> comp;
	int inComp[NodeN];


	for (int k = 0; k < NodeN; k++) {
		GSCC_vector[k] = NodeN;
	}

	for (int i = 0; i < NodeN; i++) {
		visited[i] = 0;
		inComp[i] = 0;
	}
	
	for (int i = 0; i < NodeN; i++) {
		lowlink[i] = NodeN + 1;
	}

	for (int i = 0; i < NodeN; i++) {
		if (visited[i] == 0 && del[i] == 0) {
			StrongConnect(Start_o, Adjacency, visited, del, lowlink, comp, GSCC_vector, i, index, 0, gscc_length, NodeN, tempGSCC, inComp);
		}
	}
	gscc_length++;
	if(Debugging==true){
	cout << "gscc_length  " << gscc_length << endl;
	cout << "Here is GSCC:" << endl;
	
		
		for (int k = 0; k < gscc_length; k++) {
			cout << GSCC_vector[k] << endl;		
		} 
	}
	
	if (gscc_length == 1) {
		GSCC_vector[0] = 0;
	}

	return gscc_length;
}

//*****************************************************
//--------UTILITIES----------------
/*Build the probability distribution from the final DegreeList.*/
int EmpiricalPk(int DegreeList[], int DegreeTotal,double NetP_K[], int Kmax,int NodeN){
	//First, Find the Kmax on this network:
	int HighestDegree = 0;
	
	for(int i=0;i < Kmax+1;i++){
		NetP_K[i]=0;
	}	
	
	//Count instances of each k to build P(k).
	for(int i= 0;i < NodeN;i++){
		NetP_K[DegreeList[i]]++;
		
		if(DegreeList[i] > HighestDegree) HighestDegree = DegreeList[i];
	}
		if(Debugging==true) 
			cout<<"Highest Degree in Net: "<< HighestDegree<<endl;
	
	for(int i=0;i<Kmax+1;i++){
              NetP_K[i]/=((double)NodeN);
      }
		
	return 0;
}	


void printArray(int* array, int arraylen){
	
	for(int i = 0;i< arraylen;i++){
		cout<<array[i]<<" ; ";
		if(i % 16 == 0 && i != 0){
			cout << endl;
		}
	}
	cout <<" end"<< endl;
	
	}
	
void printArray(double* array, int arraylen){
	
	for(int i = 0;i< arraylen;i++){
		cout<<array[i]<<" ; ";
		if(i % 16 == 0 && i != 0){
			cout << endl;
		}
	}
	cout <<" end"<< endl;
	
	}
	
/*Using the given arrays, assumed to be of equal length 'len', calculates the pearson correlation coefficient
 * by calculating average and square average values of each array, and using those values to generate the
 * standard deviation, and the coefficient, Rho:
 * Rho = (<xy> - <x><y>)/(sqrt(<x^2>-<x>^2)*sqrt(<y^2>-<y>^2)) * 
 * NOTE: For <y^2>-<y>^2==0 we let the correlation equal ZERO.
 * */
double PearsonCoefficient(int X[], int Y[], int len){
	double rho = 0.0;
	double covariance = 0;
	double avx = 0; double sqavx = 0;//<x>, and <x^2>
	double avy = 0; double sqavy = 0;//<y>, and <y^2>
	double sigX = 0; double sigY = 0;
	double prodav = 0;//<xy>
	for (int i = 0; i < len;i++){
		avx += X[i]; sqavx += X[i] * X[i];
		avy += Y[i]; sqavy += Y[i] * Y[i];
		prodav += X[i] * Y[i];	}
		
	avx /= (double)len; sqavx /= (double)len;
	avy /= (double)len; sqavy /= (double)len;
	prodav /= (double)len;
	
	covariance = prodav - avx*avy;
	
	if(sqavy - avy*avy==0) return 0;
	sigX = sqrt(sqavx - avx*avx);sigY = sqrt(sqavy - avy*avy);
	
	rho = covariance / (sigX * sigY);
	
	return rho;
}	
/*Using the given arrays, assumed to be of equal length 'len', calculates the pearson correlation coefficient
 * by calculating average and square average values of each array, and using those values to generate the
 * standard deviation, and the coefficient, Rho:
 * Rho = (<xy> - <x><y>)/(sqrt(<x^2>-<x>^2)*sqrt(<y^2>-<y>^2)) * 
 * NOTE: For <y^2>-<y>^2==0 we let the correlation equal ZERO.
 * */
double PearsonCoefficient(double X[], double Y[], int len){
	double rho = 0.0;
	double covariance = 0;
	double avx = 0; double sqavx = 0;//<x>, and <x^2>
	double avy = 0; double sqavy = 0;//<y>, and <y^2>
	double sigX = 0; double sigY = 0;
	double prodav = 0;//<xy>
	for (int i = 0; i < len;i++){
		avx += X[i]; sqavx += X[i] * X[i];
		avy += Y[i]; sqavy += Y[i] * Y[i];
		prodav += X[i] * Y[i];	}
		
	avx /= (double)len; sqavx /= (double)len;
	avy /= (double)len; sqavy /= (double)len;
	prodav /= (double)len;
	
	covariance = prodav - avx*avy;
	
	if(sqavy - avy*avy==0) return 0;
	sigX = sqrt(sqavx - avx*avx);sigY = sqrt(sqavy - avy*avy);
	
	rho = covariance / (sigX * sigY);
	
	return rho;
}

//*****************************************************
//---------NOISE ADDITION---------------------

//Deleting Links

void DeletingLinksUndirected (int StartTrue[], int StartFalse [], int AdjacencyTrue[], int AdjacencyFalse[], int Perm2[], int UnshuffledAdjacency[], double delta, int DegreeTotal, int NodeN) {
	int ToBeDel_1, ToBeDel_2;
	int PermIndex = 0;
	double EdgesToBeRounded = DegreeTotal * delta * 0.5;
	int EdgesToBeDel = (int) (EdgesToBeRounded + 0.5);

	for (int i = 0; i < DegreeTotal; i++) {
		AdjacencyFalse[i] = AdjacencyTrue[i];
	}

	for (int i = 0; i <= NodeN; i++) {
		StartFalse[i] = StartTrue[i];
	}

	InitAdjacencyList(StartTrue, NodeN, DegreeTotal, UnshuffledAdjacency);

	for (int i = 0; i < DegreeTotal; i++) {
		if (UnshuffledAdjacency[i] < AdjacencyFalse[i]) {
			Perm2[PermIndex] = i;
			PermIndex++;
		}
	}
	
	ShufflePerm2List(Perm2, DegreeTotal/2, EdgesToBeDel);

	for (int i = 0; i < EdgesToBeDel; i++) {
		ToBeDel_1 = AdjacencyFalse[Perm2[i]];
		AdjacencyFalse[Perm2[i]] = -1;
		ToBeDel_2 = UnshuffledAdjacency[Perm2[i]];
		for (int j = StartTrue[ToBeDel_1]; j < StartTrue[ToBeDel_1 + 1]; j++) {
			if (AdjacencyFalse[j] == ToBeDel_2) {
				AdjacencyFalse[j] = -1;
				break;
			}
		}
	}
}


void DeletingLinksDirected (int StartTrue[], int StartFalse [], int AdjacencyTrue[], int AdjacencyFalse[], int Perm2[], int DegreeList_i_False[], double delta, int DegreeTotal, int NodeN) {
	double EdgesToBeRounded = DegreeTotal * delta;
	int EdgesToBeDel = (int) (EdgesToBeRounded + 0.5);

	for (int i = 0; i < DegreeTotal; i++) {
		AdjacencyFalse[i] = AdjacencyTrue[i];
	}

	for (int i = 0; i <= NodeN; i++) {
		StartFalse[i] = StartTrue[i];
	}

	for (int i = 0; i < DegreeTotal; i++) {
		Perm2[i] = i;
	}


	ShufflePerm2List(Perm2, DegreeTotal, EdgesToBeDel);

	for (int i = 0; i < EdgesToBeDel; i++) {
		AdjacencyFalse[Perm2[i]] = -1;
		DegreeList_i_False[AdjacencyTrue[Perm2[i]]] = DegreeList_i_False[AdjacencyTrue[Perm2[i]]] - 1;
	}
}
	




int PruneForFalse(int StartFalse[], int AdjacencyFalse[], int NodeN, int& DegreeTotalFalse, int DegreeListFalse[]){
	int BadCount2 = 0;
	int oldBadCount2 = 0;
	bool jthIsBad = false;

//for i = 1:length(S)-1
	for(int i = 0;i <= NodeN-1;i++){
		oldBadCount2 = BadCount2;
	
	// for j = S(i):S(i+1)-1
		for(int j = StartFalse[i];j < StartFalse[i+1];j++){
		// Over each node's edge list.
			jthIsBad = false;
			if(AdjacencyFalse[j] == -1){
				BadCount2++;
				jthIsBad = true;
			}

			if(jthIsBad==false){
				AdjacencyFalse[j-BadCount2] = AdjacencyFalse[j];	
			} 

		}
		StartFalse[i] = StartFalse[i] - oldBadCount2;	
	}
	StartFalse[NodeN] = StartFalse[NodeN] - BadCount2;

	for (int j = DegreeTotalFalse-(BadCount2);j < DegreeTotalFalse;j++){
		AdjacencyFalse[j] = -1;
	}



	for (int i = 0; i < NodeN; i++) {
		DegreeListFalse[i] = (StartFalse[i+1] - StartFalse[i]);
	}


	DegreeTotalFalse = DegreeTotalFalse - BadCount2;


	return BadCount2;

}


// Construct EdgesToAdd List

void AddingEdgesModel1 (int Start[], int Adjacency[], int EdgesToAdd_1[], int EdgesToAdd_2[], double alpha, int NodeN, int EdgesToBeAdded, char BuildMode) {


	int rand_i;
	int rand_j;
	int EdgesAdded = 0;
	
	while (EdgesAdded < EdgesToBeAdded) {

		ChooseAgain:
		rand_i = rand() % NodeN;

		do {
			rand_j = rand() % NodeN;

			for (int i = Start[rand_i]; i < Start[rand_i + 1]; i++) {
				if (Adjacency[i] == rand_j) {
					goto ChooseAgain;
				}
			}

			for (int i = 0; i < EdgesAdded; i++) {
				if (EdgesToAdd_1[i] == rand_i && EdgesToAdd_2[i] == rand_j) {
					goto ChooseAgain;
				}
			}
		} while (rand_j == rand_i);

		EdgesToAdd_1[EdgesAdded] = rand_i;
		EdgesToAdd_2[EdgesAdded] = rand_j;

		if (BuildMode == 'u') {

			EdgesToAdd_1[EdgesAdded+1] = rand_j;
			EdgesToAdd_2[EdgesAdded+1] = rand_i;

			EdgesAdded = EdgesAdded + 2;
		}

		if (BuildMode == 'd') {
			EdgesAdded++;
		}
	}
}



void AddingEdgesModel2_Undirected (int Start[], int Adjacency[], int UnshuffledAdjacency[], int EdgesToAdd_1[], int EdgesToAdd_2[], double alpha, int DegreeTotal, int EdgesToBeAdded) {
	
	int rand_i;
	int rand_j;
	int EdgesAdded = 0;

	while (EdgesAdded < EdgesToBeAdded) {
		ChooseAgain:

		rand_i = rand() % DegreeTotal;

		do {
			rand_j = rand() % DegreeTotal;

			for (int i = Start[UnshuffledAdjacency[rand_i]]; i < Start[UnshuffledAdjacency[rand_i] + 1]; i++) {
				if (Adjacency[i] == UnshuffledAdjacency[rand_j]) {
					goto ChooseAgain;
				}
			}

			for (int i = 0; i < EdgesAdded; i++) {
				if (EdgesToAdd_1[i] == UnshuffledAdjacency[rand_i] && EdgesToAdd_2[i] == UnshuffledAdjacency[rand_j]) {
					goto ChooseAgain;
				}
			}
		} while (UnshuffledAdjacency[rand_j] == UnshuffledAdjacency[rand_i]);

		EdgesToAdd_1[EdgesAdded] = UnshuffledAdjacency[rand_i];
		EdgesToAdd_2[EdgesAdded] = UnshuffledAdjacency[rand_j];

		EdgesToAdd_1[EdgesAdded+1] = UnshuffledAdjacency[rand_j];
		EdgesToAdd_2[EdgesAdded+1] = UnshuffledAdjacency[rand_i];

		EdgesAdded = EdgesAdded + 2;
	}
}


void AddingEdgesModel2_Directed (int Start[], int Adjacency[], int UnshuffledAdjacency_i[], int UnshuffledAdjacency_o[], int EdgesToAdd_1[], int EdgesToAdd_2[], double alpha, int DegreeTotal, int EdgesToBeAdded) {
	
	int rand_i;
	int rand_j;
	int EdgesAdded = 0;

	while (EdgesAdded < EdgesToBeAdded) {
		ChooseAgain:

		rand_i = rand() % DegreeTotal;

		do {
			rand_j = rand() % DegreeTotal;

			for (int i = Start[UnshuffledAdjacency_o[rand_i]]; i < Start[UnshuffledAdjacency_o[rand_i] + 1]; i++) {
				if (Adjacency[i] == UnshuffledAdjacency_i[rand_j]) {
					goto ChooseAgain;
				}
			}

			for (int i = 0; i < EdgesAdded; i++) {
				if (EdgesToAdd_1[i] == UnshuffledAdjacency_o[rand_i] && EdgesToAdd_2[i] == UnshuffledAdjacency_i[rand_j]) {
					goto ChooseAgain;
				}
			}
		} while (UnshuffledAdjacency_i[rand_j] == UnshuffledAdjacency_o[rand_i]);

		EdgesToAdd_1[EdgesAdded] = UnshuffledAdjacency_o[rand_i];
		EdgesToAdd_2[EdgesAdded] = UnshuffledAdjacency_i[rand_j];
		EdgesAdded++;
	}
}
	
//QuickSort Code
// If this code is used to sort only one array, the second input array (EdgesToAdd_2) will not be changed
int Partition_1(int EdgesToAdd_1[], int EdgesToAdd_2[], int p, int r)

{
    int pivot = EdgesToAdd_1[r];

    while ( p < r )
    {
        while ( EdgesToAdd_1[p] < pivot )
            p++;

        while ( EdgesToAdd_1[r] > pivot )
            r--;

        if ( EdgesToAdd_1[p] == EdgesToAdd_1[r] )
            p++;
        else if ( p < r )
        {
            int tmp_1 = EdgesToAdd_1[p];
			int tmp_2 = EdgesToAdd_2[p];

			EdgesToAdd_1[p] = EdgesToAdd_1[r];
			EdgesToAdd_1[r] = tmp_1;

			EdgesToAdd_2[p] = EdgesToAdd_2[r];
			EdgesToAdd_2[r] = tmp_2;

        }
    }

    return r;
}





void Quicksort_1(int EdgesToAdd_1[], int EdgesToAdd_2[], int p, int r)
{
    if ( p < r )
    {
        int j = Partition_1(EdgesToAdd_1, EdgesToAdd_2, p, r);        
        Quicksort_1(EdgesToAdd_1, EdgesToAdd_2, p, j-1);
        Quicksort_1(EdgesToAdd_1, EdgesToAdd_2, j+1, r);

    }

}



int Partition_2(double KeyValues[], int SortOrder[], int p, int r)
{
    double pivot = KeyValues[r];

    while ( p < r )
    {
        while ( KeyValues[p] < pivot )
            p++;

        while ( KeyValues[r] > pivot )
            r--;

        if ( KeyValues[p] == KeyValues[r] )
            p++;
        else if ( p < r )
        {
            double tmp_1 = KeyValues[p];
			int tmp_2 = SortOrder[p];

			KeyValues[p] = KeyValues[r];
			KeyValues[r] = tmp_1;

			SortOrder[p] = SortOrder[r];
			SortOrder[r] = tmp_2;
        }
    }

    return r;
}


void Quicksort_2(double KeyValues[], int SortOrder[], int p, int r)
{
    if ( p < r )
    {
        int j = Partition_2(KeyValues, SortOrder, p, r);        
        Quicksort_2(KeyValues, SortOrder, p, j-1);
        Quicksort_2(KeyValues, SortOrder, j+1, r);
    }
}


// Add Edges From EdgesToAdd List

void AddEdgesToStartAndAdjacency_Undirected (int FinalAdjacencyFalse[], int FinalStartFalse[], int AdjacencyFalse[], int StartFalse[], int DegreeListFalse[], int EdgesToBeAdded, int& DegreeTotalFalse, int NodeN, int EdgesToAdd_1[], int EdgesToAdd_2[]) {

	int index = 0;
	int old_m = 0;
	int new_m = 0;

	//Assign FinalStartFalse
		for (int i = 0; i <= NodeN; i++) {
			FinalStartFalse[i] = StartFalse[i];
		}


	for (int k = 0; k < EdgesToBeAdded; k++) {
		new_m = EdgesToAdd_1[k];

		
		for (int i = old_m; i < new_m; i++) {
			
			for (int j = StartFalse[i]; j < StartFalse[i + 1]; j++) {
				FinalAdjacencyFalse[j + index] = AdjacencyFalse[j];
			}
		}
		
		FinalAdjacencyFalse[StartFalse[new_m]+index] = EdgesToAdd_2[k];
		index++;
		old_m = new_m;
	
		for (int z = new_m+1; z <= NodeN; z++) {
			FinalStartFalse[z]++;
		}
	}

	for (int j = StartFalse[new_m]; j < DegreeTotalFalse; j++) {
		FinalAdjacencyFalse[j + index] = AdjacencyFalse[j];
	}

	DegreeTotalFalse += EdgesToBeAdded;

	for (int i = 0; i < NodeN; i++) {
		DegreeListFalse[i] = FinalStartFalse[i + 1] - FinalStartFalse[i];
	}
}



void AddEdgesToStartAndAdjacency_Directed (int FinalAdjacencyFalse[], int FinalStartFalse_o[], int FinalStartFalse_i[], int AdjacencyFalse[], int StartFalse_o[], int DegreeList_i_False[], int DegreeList_o_False[], int EdgesToBeAdded, int& DegreeTotalFalse, int NodeN, int EdgesToAdd_1[], int EdgesToAdd_2[]) {

	int index = 0;
	int old_m = 0;
	int new_m = 0;

	//Assign FinalStartFalse
		for (int i = 0; i <= NodeN; i++) {
			FinalStartFalse_o[i] = StartFalse_o[i];
		}


	for (int k = 0; k < EdgesToBeAdded; k++) {
		new_m = EdgesToAdd_1[k];

		
		for (int i = old_m; i < new_m; i++) {
			
			for (int j = StartFalse_o[i]; j < StartFalse_o[i + 1]; j++) {
				FinalAdjacencyFalse[j + index] = AdjacencyFalse[j];
			}
		}
		
		FinalAdjacencyFalse[StartFalse_o[new_m]+index] = EdgesToAdd_2[k];
		index++;
		old_m = new_m;
	
		for (int z = new_m+1; z <= NodeN; z++) {
			FinalStartFalse_o[z]++;
		}
	}
	for (int j = StartFalse_o[new_m]; j < DegreeTotalFalse; j++) {
		FinalAdjacencyFalse[j + index] = AdjacencyFalse[j];
	}

	DegreeTotalFalse += EdgesToBeAdded;

	for (int i = 0; i < NodeN; i++) {
		DegreeList_o_False[i] = FinalStartFalse_o[i + 1] - FinalStartFalse_o[i];
	}

	for (int i = 0; i < EdgesToBeAdded; i++) {
		DegreeList_i_False[EdgesToAdd_2[i]]++;
	}
	
	FinalStartFalse_i[0] = 0;
	for (int i = 1; i <= NodeN; i++) {
		FinalStartFalse_i[i] = DegreeList_i_False[i - 1] + FinalStartFalse_i[i-1];
	}
		
}


//-----------ATTACK ALGORITHMS--------------
/*Rank all of the nodes in this network according to degree, using a counting
 * sort.(pseudocode per wikipedia.)*/
int DegreeCentralityRanking(int DegreeList[], int NodeN, int DegCentrality[]){
	int count[NodeN];
	int DegCentrality1[NodeN];
	int total = 0; int oldcount = 0;
	
	for(int i=0;i<NodeN;i++) count[i] = 0;
	
	for(int i=0;i<NodeN;i++) count[DegreeList[i]]++;//count occurances of nodes of each deg.
		
	for(int i=0;i<NodeN;i++){
		//turn the count list into an index list.
		oldcount = count[i];
		count[i] = total;
		total += oldcount;
		}
		
	for(int i=0;i<NodeN;i++)	{
		DegCentrality1[count[DegreeList[i]]] = i;
		count[DegreeList[i]]++;
		}
		
	for(int i=0;i<NodeN;i++) DegCentrality[i] = DegCentrality1[NodeN-i-1];	
		
	if(Debugging==true){cout<<"Ordered Degree list(descending):"<<endl; printArray(DegCentrality,NodeN);}	
	
	return 0;
}

void QuicksortCentralityRanking (double KeyValues[], int NodeN, int SortOrder[]) {
	for (int i = 0; i < NodeN; i++) {
		SortOrder[i] = i;
	}
	Quicksort_2(KeyValues, SortOrder, 0, NodeN-1);
}


/*Using list of nodes deleted so far, and the Intact Start/Adjacency/Degreelist, a new
 the DegreeList is recalculated for re-ranking of nodes.  */
int DegreeListUpdate(int OrigDegreeList[], int Start[], int NodeN, int Adjacency[], int DegreeTotal, int del[], int NewDegList[])
{
	for(int i = 0;i<NodeN;i++){NewDegList[i] = OrigDegreeList[i];}
	//copy the original list.
	
//	NewDegList[DeletedIndex]=0;
//	for(int Ai = Start[DeletedIndex];Ai<Start[DeletedIndex+1];Ai++)
//			NewDegList[Adjacency[Ai]]--;

	for(int i = 0;i<NodeN;i++){
		if(del[i]==1){
		NewDegList[i] = 0;
		//for deleted nodes, clear their degcount. 
		//decrement the degcount for all connected nodes by one. 
			for(int Ai = Start[i];Ai<Start[i+1];Ai++)
				if(del[Adjacency[Ai]]==0) NewDegList[Adjacency[Ai]]--;
		}
	}
	//printArray(NewDegList,NodeN);	
	return 0;
}

/*FOR NETS WITH DELETED NODES, Run the same Centrality contribution algorithm as kept above only skip BFS'ing 
 * out to any node 'i' for which del[i] == 1.*/
int BetCentContrib(int InStart[],int OutStart[], int NodeN, int Adjacency[], int DegreeTotal, double Bpart[], int root, int del[]){

// q: a queue of nodes to visit
// qs: start of used portion of q
// qe: end of used portion of q
	int Q[NodeN];
	int qs=0;int qe=0;
	
//distance to a node
	int DtoRoot[NodeN];
	
//number of paths of minimum distance to each node
	int PathsFromRoot[NodeN];
  
//!Initialize queue to contain only root node, distances to be max 
//!possible, and betweennesses to be 0.
 for(int n = 0;n < NodeN;n++){
	 Q[n]=NodeN+1;
	 DtoRoot[n] = NodeN+1;
	 Bpart[n] = 0;
	 PathsFromRoot[n]=0;
 }
 Q[0]=root;
 DtoRoot[root] = 0;
 PathsFromRoot[root] = 1; 
 
int parents[DegreeTotal];
for(int p = 0; p< DegreeTotal;p++){
	parents[p] = NodeN+1;}

//A Breadth First Search for the shortest paths to all other nodes that begin with 
//the root given. 
 int i, j, IndIncr;
 do 
 {
	i = Q[qs];
	qs = qs + 1;
		//Directed Version; branch out uses InStart, parent check uses OutStart
		for(int e = OutStart[i];e<OutStart[i+1];e++){	
			j = Adjacency[e];
			if(del[j]==1) continue;//only BFS out toward nodes that are present.
			if(DtoRoot[j] >= (DtoRoot[i]+1)){//condition for a shortest path, possibly non-unique
				if(DtoRoot[j] > DtoRoot[i]+1){//condition for unexplored node
					qe=qe+1;
					Q[qe] = j;
					DtoRoot[j] = DtoRoot[i] + 1;
				}
				PathsFromRoot[j] = PathsFromRoot[j] + PathsFromRoot[i];
				
			//For that node's particular parent list, only add if no shortest path
			//parent occupies that position already.
					IndIncr = 0;
					while(parents[InStart[j]+IndIncr] < NodeN+1 ){IndIncr++; }
					parents[InStart[j]+IndIncr] = i;
				}			
		}
//cout<<"Paths from root:"<<endl;
//printArray(PathsFromRoot, NodeN);		
//cout<<"Distances to r("<<root<<"):"<<endl;
//printArray(DtoRoot, NodeN);	
//cout<<"Queue("<<qs<<"..."<<qe<<"):"<<endl;
//printArray(Q, NodeN);
//cout<<"Parents:"<<endl;
//printArray(parents,DegreeTotal);
//cout<<"-->Next:"<<Q[qs]<<endl;
	 
}while (qs <= qe);

//------------------- 

//Betweenness, given num paths to root, distances to root, and paths via parent trace-back;
//calculates Betweenness Centrality Contribution for all nodes != to the root.
int k;
for(int Queuei = qe; Queuei>= 0; Queuei--){
	IndIncr = 0;

	k=Q[Queuei];

		while(parents[InStart[k]+IndIncr]<NodeN+1 && InStart[k]+IndIncr < InStart[k+1]){
	//pull parents from the part of parents list that is filled, but not going over into the
	//next node's parents.
		j = parents[InStart[k]+IndIncr];
		Bpart[j] = Bpart[j]+ (1.0 + Bpart[k] )*(double)PathsFromRoot[j] / (double)PathsFromRoot[k];
		IndIncr++;
	}		
}  
 Bpart[root] = 0; 
 
 if(Debugging==true)
	printArray(Bpart,NodeN);
 //working backward from the end of the queue. Algorithm per S.S.	
//INITIALIZE: bb(i)=0 for all i
//for each i (going backwards in the queue)
//for each parent of i,
//       bb(parent)=bb(parent)+(1.+bb(i))*((float)np(parent))/((float)np(i))
//AFTERWARDS: bb(root node)=0.
 return 0;

}

/*For both directed and undirected networks, calculates the betweenness of each node by counting the number of times
 * each nodes occurs in the shortest paths between all pairs(weighting when there is more than 1 s.p.). 
 * IN: In and Out Start lists for Directed. FOR UNDIRECTED, SEND 'START' for both inputs. 
 * 'Btotal' will have Betweenness values upon return.  * */
 double *Betweenness(int InStart[], int OutStart[], int NodeN, int Adjacency[], int DegreeTotal, double Btotal[],int del[]){
	double Bpart[NodeN];
	//every Node will run BetCentContrib, using shortest paths from that node to all
	//other nodes to determine the contribution of paths from that node to the betweenness
	//of every other node.	
	for(int i=0;i<NodeN;i++) Btotal[i]=0.;
	for(int root = 0;root < NodeN;root++)
	{
		
		if(del[root] == 1) continue;//only use undeleted nodes for a root.
		
		for(int i = 0;i< NodeN;i++){
		Bpart[i] = 0.0;
		}				
		
		if(Debugging==true)	cout<<"Contrib, Node "<<root<<endl;
		BetCentContrib(InStart, OutStart, NodeN, Adjacency, DegreeTotal, Bpart, root, del);

		for(int j = 0;j<NodeN;j++){
		Btotal[j] = Btotal[j] + Bpart[j];
		}
		
		if(root % 1000 == 0){
				cout<<"* ";
					if(root % 15000==0){cout<<endl;}
						}
	}
	 
	if(Debugging==true)
		{cout<<"Betweenness:"<<endl; printArray(Btotal, NodeN);}
	
	return Btotal;
}

/*VS 'centrality' measure version of these algorithms, only the 'del' check should be different.
 * From the Start and Adjacency lists, determine the eigenvector corresponding to the largest
 * eigenvalue of this matrix. The 'eigVector' passed will have the vector filled upon return. 
The largest eigenvalue is also returned as a double by the function. 
 * FOR DIRECTED NETWORKS: Pass-in the out_Start list.(From the derivation of the Dynamical 
 * Importance approx., Restrepo,et al 2006: "We consider a network as a directed graph with N nodes,
and we associate to it a N x N adjacency matrix whose
elements Aij are positive if there is a link going from node i
to node j with i != j and zero otherwise (Aii = 0)."		 */
double eigenVector(int Start[], int NodeN, int Alist[], int DegreeTotal,double eigVector[],bool leftVect, int del[]){
	
	//Prepare a pair of vectors of length N, for iteration. Let the iteration begin with 
	//a random vector, which is then normalized.
	double newV[NodeN];	
	double oldV[NodeN];
	double NormSum = 0.0;
	
	int a = 0;
	int power = 0;	
	
	double nsL=0;//New and Old Shifted Lambdas(\lambda + 1 = nsL)
	double osL=0;
	double Lambda = 0;
	
	
	//init a blank V_n+1, let V_0 be a random vector, which is normalized.
	for(int i = 0;i<NodeN;i++) {
		newV[i] = 0.0;
		if(del[i]==0) oldV[i] = Pdist(rng);
		else oldV[i] = 0.;
		NormSum += oldV[i]*oldV[i];
	}
		
	NormSum = sqrt(NormSum);		
	for(int i = 0;i < NodeN;i++) oldV[i] /= NormSum;	
	osL = NormSum;
	
	if(Debugging==true){
	cout<<"V_"<< 0 << endl;
	printArray(oldV,NodeN);	}

do{//while lambda has not converged to within 1e-5.
	
	power++;
	//Multiply the Adjacency matrix to oldV to give newV.
	for(int i =0; i < NodeN;i++){//V_n+1 = V_n + A*V_n
		newV[i] = oldV[i];
	}
	for(int i=0; i<NodeN; i++){
		if(del[i]==1) {continue;}
		for(int j=Start[i];j<Start[i+1];j++){
			a = Alist[j];
			if (del[a]==0){
				if(leftVect==false)newV[i] += oldV[a];
				else newV[a] += oldV[i];
			}
		//Per notes: The left eigenvector of A is equal to the right 
		//	eigenvector of A^T(transpose). The switch of indices uses the transpose.  
		}	
	}	
	
	//Calculate the Norm(root(Sum(V)) for n+1.
	for(int i =0; i<NodeN;i++) nsL += newV[i]*newV[i];
	nsL = sqrt(nsL);
	
	//normalize the new V_n+1
	for(int i =0; i < NodeN;i++) oldV[i] = newV[i] / nsL;
	//cout << power << "   " << nsL << "   " << newV[0] << "   " << nsL-osL << endl;
	
	//cout<<"nsL = " << nsL <<". V_"<< power <<": "<< endl;
	//printArray(oldV,NodeN);	
	
	//Check for convergence of the eigenvalue; The magnitude of the V_n+1 we obtain.	
	if(fabs(nsL - osL) < 1e-5) break;
	else {	osL = nsL;
			nsL = 0;		}
	
}while(true);

	Lambda = nsL - 1.0;	
	if(Debugging==true)	cout<<"Found EigenValue- "<<Lambda<<endl;
	
	if(Debugging==true){ cout<<"Found EigenVector: \n";printArray(oldV,NodeN);}
	for(int p =0;p<NodeN;p++) eigVector[p] = oldV[p];	
	
	return Lambda;
	
}

/*Using the approximation given in Restrepo, et al. 2006:
 * I(k) = u_k*v_k/(u^T * v ) The DI of the kth node is approximated by taking the product
 * of the kth elements of the left and right eigenVectors and normalizing them by the dot product
 * of those vectors. The given array, DynImp[] will be filled with these values.  */
int DynamicalImportance_Directed(int Start[], int Adjacency[], int DegreeTotal, int NodeN, double DynImp[],int del[]){
	//Take the dot product of the left/right eigenvectors. This is the normalization factor.
	
	double lVect[NodeN];
	double rVect[NodeN];
	eigenVector(Start, NodeN, Adjacency, DegreeTotal, lVect, true, del);
	eigenVector(Start, NodeN, Adjacency, DegreeTotal, rVect, false, del);
	
	double NormFactor= 0.0;	
	for(int i = 0;i<NodeN;i++){	NormFactor += lVect[i] * rVect[i];	}
	//Note: Should be 1.0 for the undirected case; as the eigVector should be normalized before call.
	
	//Into DynImp, calculate the approximate DI for every node, as given above.
	for(int i = 0;i<NodeN;i++){ DynImp[i] = lVect[i] * rVect[i] / NormFactor;	}	
	
	return 0;
}


int DynamicalImportance_Undirected(int Start[], int Adjacency[], int DegreeTotal, int NodeN, double DynImp[],int del[]){
	
	double eVect[NodeN];	
	eigenVector(Start, NodeN, Adjacency, DegreeTotal, eVect, false, del);
	
	double NormFactor= 0.0;	
	for(int i = 0;i<NodeN;i++){	
		NormFactor += eVect[i] * eVect[i];	
		}
	
	//Into DynImp, calculate the approximate DI for every node, as given above.
	for(int i = 0;i<NodeN;i++){
	DynImp[i] = eVect[i] * eVect[i] / NormFactor;
	}	

	return 0;
}


//*****************************************************
