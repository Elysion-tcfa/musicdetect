#include <cstdio>
#include<iostream>
#include "RTree.h"
#include<time.h>
#include<string>
#include<fstream>
#include<tchar.h>
#include "config.h"


using namespace std;
#define K 10
#define T 25
struct Rect
{
  Rect()  {}

  


  double min[DIM];
  double max[DIM];
};

int localrecord[450]={0};
int hitcount=0;
int havebuild=0;
bool MySearchCallback(int id, void* arg) 
{
	localrecord[hitcount++]=id;
  return true; // keep going
}

double * paa(double* original, int length)
{
	double * now = new double[DIM];  //开辟数组空间，长度为M
	
	for(int i = 1;i <= DIM;i++)   
	{
		double num = 0;
		int j = length * (i - 1) / DIM+ 1;
		int max = length * i / DIM;
		for(;j <= max;j++)
		{
			num += original[j-1];
		}
		now[i-1] = DIM * num / length;
		//cout << now[i-1] << endl;
	}
	return now;
}
RTree<int, double, DIM, double> tree;//type of id, type of coordinate, number of dim, type of volume result




extern "C" _declspec(dllexport) int getsearch(int record[],char * filename )
{
	
	fstream file;
	file.open(filename, ios::in);
	RTFileStream rts;//creat a file stream with necessary check, or tree.Load directly is also ok.
	hitcount=0;
	if(havebuild==0)
		tree.Load("h:\\projectdata\\database3.rtr");
	havebuild=1;
	double value;
	int length;
	Rect searchrec;
	 file>>length;
  for(int j=0;j<length;j++)
	  file>>value;
  file>>length;
	  for(int j=0;j<length;j++)
	  {
		  
		 file>>value;
		  searchrec.min[j]=value;
		  searchrec.max[j]=value;
	  }
	  for(int i=0;i<length;i++)
	  {
		  int tempmax=-99,tempmin=99;
		  for(int j=max(i-K,0);j<=min(length-1,i+K);j++)
		  {
			  if(searchrec.min[i]<tempmin)
				  tempmin=searchrec.min[i];
			  if(searchrec.max[i]>tempmax)
				  tempmax=searchrec.max[i];
		  }
		  searchrec.max[i]=tempmax;
		  searchrec.min[i]=tempmin;
	  }
	    file.close();
 // double * searchseq=paa(des,length);
  /*Rect searchrec;
  for(int i=0;i<DIM;i++)
  {
	  searchrec.max[i]=T;
	  searchrec.min[i]=-T;
  }*/
		for(int j=0;j<length;j++)
	  {
		  
		
		  searchrec.min[j]-=T;
		  searchrec.max[j]+=T;
	  }
  tree.Search(searchrec.min,searchrec.max,MySearchCallback,NULL);
  for(int i=0;i<hitcount;i++)
	  record[i]=localrecord[i];

  return hitcount;




 
}