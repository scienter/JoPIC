#include <stdio.h>
#include <stdlib.h>
#include "mesh.h"
#include "mpi.h"

void MPI_TransferF_Zminus(Domain *D,int numField)
{
  int n,i,j,k,numberData,start,end,num,nx,ny,share;
  int istart,iend,jstart,jend,kstart,kend;
  int myrank, nTasks,rank; 
  double *data;

  MPI_Status status;         

  istart=D->istart;    iend=D->iend;
  jstart=D->jstart;    jend=D->jend;
  kstart=D->kstart;    kend=D->kend;
  nx=1;	ny=1;
  if(D->dimension>1) { nx=iend+3; } else ;
  if(D->dimension>2) { ny=jend+3; } else ;

  MPI_Comm_size(MPI_COMM_WORLD, &nTasks);
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);            

  rank=(myrank%(D->M*D->N))/D->M;

  share=D->shareZ;
  num=numField*nx*ny*share;
  data = (double *)malloc(num*sizeof(double ));

  //Transferring even ~ odd cores 
  start=0; 
  for(n=0; n<numField; n++)
    for(i=0; i<nx; i++)
      for(j=0; j<ny; j++)	
      {
        for(k=0; k<share; k++) 
          data[start+k]=D->shareF[n][i][j][k+kstart]; 
        start+=share; 
      }

  if(rank%2==0 && rank!=D->N-1)
  {
    MPI_Recv(data,num,MPI_DOUBLE,D->nextZrank,D->nextZrank, MPI_COMM_WORLD,&status);  
    start=0; 
    for(n=0; n<numField; n++)
      for(i=0; i<nx; i++)
        for(j=0; j<ny; j++)	
        {
          for(k=0; k<share; k++) 
            D->shareF[n][i][j][k+kend]=data[start+k]; 
          start+=share; 
        }
  }
  else if(rank%2==1)
    MPI_Send(data,num,MPI_DOUBLE,D->prevZrank,myrank,MPI_COMM_WORLD);             
  MPI_Barrier(MPI_COMM_WORLD);

  //Transferring odd ~ even cores 
  start=0; 
  for(n=0; n<numField; n++)
    for(i=0; i<nx; i++)
      for(j=0; j<ny; j++)	
      {
        for(k=0; k<share; k++) 
          data[start+k]=D->shareF[n][i][j][k+kstart]; 
        start+=share; 
      }

  if(rank%2==1 && rank!=D->N-1)
  {
    MPI_Recv(data,num,MPI_DOUBLE,D->nextZrank,D->nextZrank, MPI_COMM_WORLD,&status);  
    start=0; 
    for(n=0; n<numField; n++)
      for(i=0; i<nx; i++)
        for(j=0; j<ny; j++)	
        {
          for(k=0; k<share; k++) 
            D->shareF[n][i][j][k+kend]=data[start+k]; 
          start+=share; 
        }
  }
  else if(rank%2==0 && rank!=0)
    MPI_Send(data,num,MPI_DOUBLE,D->prevZrank,myrank,MPI_COMM_WORLD);             
  MPI_Barrier(MPI_COMM_WORLD);
  free(data);
}

void MPI_TransferF_Zplus(Domain *D,int numField)
{
  int n,i,j,k,num,nx,ny,share;
  int istart,iend,jstart,jend,kstart,kend;
  int myrank, nTasks,rank,start; 
  double *data;

  MPI_Status status;         
  MPI_Comm_size(MPI_COMM_WORLD, &nTasks);
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);            
   
  istart=D->istart;    iend=D->iend;
  jstart=D->jstart;    jend=D->jend;
  kstart=D->kstart;    kend=D->kend;
  nx=1;	ny=1;
  if(D->dimension>1) { nx=iend+3; } else ;
  if(D->dimension>2) { ny=jend+3; } else ;

  share=D->shareZ;
  rank=(myrank%(D->M*D->N))/D->M;
  num=numField*nx*ny*(share-1);
  data = (double *)malloc(num*sizeof(double ));

  //Transferring even ~ odd cores 
  if(rank%2==1)
  {
    MPI_Recv(data,num,MPI_DOUBLE,D->prevZrank,D->prevZrank, MPI_COMM_WORLD,&status);  
    start=0;
    for(n=0; n<numField; n++)
      for(i=0; i<nx; i++)
        for(j=0; j<ny; j++)	
        {
          for(k=1; k<share; k++) 
            D->shareF[n][i][j][kstart-k]=data[start+k-1]; 
          start+=share-1; 
        }
  }
  else if(rank%2==0 && rank!=D->N-1)
  {
    start=0; 
    for(n=0; n<numField; n++)
      for(i=0; i<nx; i++)
        for(j=0; j<ny; j++)	
        {
          for(k=1; k<share; k++) 
            data[start+k-1]=D->shareF[n][i][j][kend-k]; 
          start+=share-1; 
        }
      
    MPI_Send(data,num,MPI_DOUBLE,D->nextZrank,myrank, MPI_COMM_WORLD);             
  }
  MPI_Barrier(MPI_COMM_WORLD);

  //Transferring odd ~even ~ cores 
  if(rank%2==0 && rank!=0)
  {
    MPI_Recv(data,num,MPI_DOUBLE,D->prevZrank,D->prevZrank, MPI_COMM_WORLD,&status);  
    start=0;
    for(n=0; n<numField; n++)
      for(i=0; i<nx; i++)
        for(j=0; j<ny; j++)	
        {
          for(k=1; k<share; k++) 
            D->shareF[n][i][j][kstart-k]=data[start+k-1]; 
          start+=share-1; 
        }
  }
  else if(rank%2==1 && rank!=D->N-1)
  {
    start=0; 
    for(n=0; n<numField; n++)
      for(i=0; i<nx; i++)
        for(j=0; j<ny; j++)	
        {
          for(k=1; k<share; k++) 
            data[start+k-1]=D->shareF[n][i][j][kend-k]; 
          start+=share-1; 
        }
      
    MPI_Send(data,num,MPI_DOUBLE,D->nextZrank,myrank, MPI_COMM_WORLD);             
  }
  MPI_Barrier(MPI_COMM_WORLD);

  free(data);
}

void MPI_TransferF_Period_Z(Domain *D,int numField)
{
  int n,i,j,k,numberData,start,end,num,nx,ny,share;
  int istart,iend,jstart,jend,kstart,kend;
  int myrank, nTasks,rank; 
  double *data;

  MPI_Status status;         

  istart=D->istart;    iend=D->iend;
  jstart=D->jstart;    jend=D->jend;
  kstart=D->kstart;    kend=D->kend;
  nx=1;	ny=1;
  if(D->dimension>1) { nx=iend+3; } else ;
  if(D->dimension>2) { ny=jend+3; } else ;

  MPI_Comm_size(MPI_COMM_WORLD, &nTasks);
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);            

  rank=(myrank%(D->M*D->N))/D->M;
  share=D->shareZ;

  // Bottom to Up 
  num=numField*nx*ny*share;
  data = (double *)malloc(num*sizeof(double ));

  start=0; 
  for(n=0; n<numField; n++)
    for(i=0; i<nx; i++)
      for(j=0; j<ny; j++)	
      {
        for(k=0; k<share; k++) 
          data[start+k]=D->shareF[n][i][j][kstart+k]; 
        start+=share; 
      }

  if(D->N==1) {
    start=0; 
    for(n=0; n<numField; n++)
      for(i=0; i<nx; i++)
        for(j=0; j<ny; j++)	
        {
          for(k=0; k<share; k++) 
            D->shareF[n][i][j][kend+k]=data[start+k]; 
          start+=share; 
        }
	} else {
	  if(rank==0)
       MPI_Send(data,num,MPI_DOUBLE,myrank+(D->N-1)*D->M,myrank,MPI_COMM_WORLD);             
     else if(rank==D->N-1) 	{
	    MPI_Recv(data,num,MPI_DOUBLE,myrank-(D->N-1)*D->M,myrank-(D->N-1)*D->M, MPI_COMM_WORLD,&status);  
       start=0; 
       for(n=0; n<numField; n++)
         for(i=0; i<nx; i++)
           for(j=0; j<ny; j++)	
           {
             for(k=0; k<share; k++) 
               D->shareF[n][i][j][kend+k]=data[start+k]; 
             start+=share; 
           }
     }
     MPI_Barrier(MPI_COMM_WORLD);
  }
  free(data);

  // Up to Bottom 
  num=numField*nx*ny*(share-1);
  data = (double *)malloc(num*sizeof(double ));

  start=0; 
  for(n=0; n<numField; n++)
    for(i=0; i<nx; i++)
      for(j=0; j<ny; j++)	
      {
        for(k=1; k<share; k++) 
          data[start+k-1]=D->shareF[n][i][j][kend-k]; 
        start+=share-1; 
      }

  if(D->N==1) {   
    start=0;
    for(n=0; n<numField; n++)
      for(i=0; i<nx; i++)
        for(j=0; j<ny; j++)	
        { 	
          for(k=1; k<share; k++) 
            D->shareF[n][i][j][kstart-k]=data[start+k-1]; 
          start+=share-1; 
        }
	} else {
	  if(rank==D->N-1)
       MPI_Send(data,num,MPI_DOUBLE,myrank-(D->N-1)*D->M,myrank, MPI_COMM_WORLD);             
     else if(rank==0)  {
       MPI_Recv(data,num,MPI_DOUBLE,myrank+(D->N-1)*D->M,myrank+(D->N-1)*D->M, MPI_COMM_WORLD,&status);  
       start=0;
       for(n=0; n<numField; n++)
         for(i=0; i<nx; i++)
           for(j=0; j<ny; j++)	
           {
             for(k=1; k<share; k++) 
               D->shareF[n][i][j][kstart-k]=data[start+k-1]; 
             start+=share-1; 
           }
     }
     MPI_Barrier(MPI_COMM_WORLD);
  }
  free(data);
}



//------------------------------ Current, Density ------------------------------------------//


void MPI_TransferJ_Zminus(Domain *D,int numField)
{
  int n,i,j,k,start,num,nx,ny,share;
  int istart,iend,jstart,jend,kstart,kend;
  int myrank, nTasks,rank; 
  double *data,tmp;

  MPI_Status status;         

  istart=D->istart;    iend=D->iend;
  jstart=D->jstart;    jend=D->jend;
  kstart=D->kstart;    kend=D->kend;
  nx=1;	ny=1;
  if(D->dimension>1) { nx=iend+3; } else ;
  if(D->dimension>2) { ny=jend+3; } else ;

  MPI_Comm_size(MPI_COMM_WORLD, &nTasks);
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);            

  share=D->shareZ;
  rank=(myrank%(D->M*D->N))/D->M;
  num=numField*nx*ny*(share-1);
  data = (double *)malloc(num*sizeof(double ));

  //Transferring even ~ odd cores 


  if(rank%2==0 && rank!=D->N-1)
  {
    MPI_Recv(data,num,MPI_DOUBLE,D->nextZrank,D->nextZrank, MPI_COMM_WORLD,&status);  
    start=0; 
    for(n=0; n<numField; n++)
      for(i=0; i<nx; i++)
        for(j=0; j<ny; j++) 
		  {	
          for(k=1; k<share; k++) 
            D->shareF[n][i][j][kend-k]+=data[start+k-1]; 
          start+=share-1; 
        }
  }
  else if(rank%2==1) 
  {
    start=0; 
    for(n=0; n<numField; n++)
      for(i=0; i<nx; i++)
        for(j=0; j<ny; j++) 
        {
          for(k=1; k<share; k++) {
            data[start+k-1]=D->shareF[n][i][j][kstart-k];
            D->shareF[n][i][j][kstart-k]=0.0;
          }
          start+=share-1; 
        }    
    MPI_Send(data,num,MPI_DOUBLE,D->prevZrank,myrank,MPI_COMM_WORLD);  
  }           
  MPI_Barrier(MPI_COMM_WORLD);

  //Transferring odd ~ even cores                     
  if(rank%2==1 && rank!=D->N-1)
  {
    MPI_Recv(data,num,MPI_DOUBLE,D->nextZrank,D->nextZrank, MPI_COMM_WORLD,&status);  
    start=0; 
    for(n=0; n<numField; n++)
      for(i=0; i<nx; i++)
        for(j=0; j<ny; j++)
		  {
          for(k=1; k<share; k++) 
            D->shareF[n][i][j][kend-k]+=data[start+k-1]; 
          start+=share-1; 
        }
  }
  else if(rank%2==0 && rank!=0) 
  {
    start=0; 
    for(n=0; n<numField; n++)
      for(i=0; i<nx; i++)
        for(j=0; j<ny; j++)
        {
          for(k=1; k<share; k++) {
            data[start+k-1]=D->shareF[n][i][j][kstart-k];
            D->shareF[n][i][j][kstart-k]=0.0;
          }
          start+=share-1; 
        }
    MPI_Send(data,num,MPI_DOUBLE,D->prevZrank,myrank, MPI_COMM_WORLD);
  }             
  MPI_Barrier(MPI_COMM_WORLD);
    
  free(data);
}

void MPI_TransferJ_Zplus(Domain *D,int numField)
{
  int n,i,j,k,start,num,nx,ny,share;
  int istart,iend,jstart,jend,kstart,kend;
  int myrank, nTasks,rank; 
  double *data,tmp;

  MPI_Status status;         
   
  istart=D->istart;    iend=D->iend;
  jstart=D->jstart;    jend=D->jend;
  kstart=D->kstart;    kend=D->kend;
  nx=1;	ny=1;
  if(D->dimension>1) { nx=iend+3; } else ;
  if(D->dimension>2) { ny=jend+3; } else ;

  MPI_Comm_size(MPI_COMM_WORLD, &nTasks);
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);            

  share=D->shareZ;
  rank=(myrank%(D->M*D->N))/D->M;
  num=numField*nx*ny*share;
  data = (double *)malloc(num*sizeof(double ));

  //Transferring even ~ odd cores       
  if(rank%2==1)
  {
    MPI_Recv(data,num,MPI_DOUBLE,D->prevZrank,D->prevZrank, MPI_COMM_WORLD,&status);  
    start=0;
    for(n=0; n<numField; n++)
      for(i=0; i<nx; i++)
        for(j=0; j<ny; j++)
        {
          for(k=0; k<share; k++) 
            D->shareF[n][i][j][kstart+k]+=data[start+k]; 
          start+=share; 
        }
  }
  else if(rank%2==0 && rank!=D->N-1) {
    start=0; 
    for(n=0; n<numField; n++)
      for(i=0; i<nx; i++)
        for(j=0; j<ny; j++)
        {
          for(k=0; k<share; k++) {
            data[start+k]=D->shareF[n][i][j][kend+k];
            D->shareF[n][i][j][kend+k]=0.0;
          }
          start+=share; 
        }
    MPI_Send(data,num,MPI_DOUBLE,D->nextZrank,myrank, MPI_COMM_WORLD);  
  }           
  MPI_Barrier(MPI_COMM_WORLD);

  //Transferring odd ~ even cores                     
  if(rank%2==0 && rank!=0)
  {
    MPI_Recv(data,num,MPI_DOUBLE,D->prevZrank,D->prevZrank,MPI_COMM_WORLD,&status);  
    start=0;
    for(n=0; n<numField; n++)
      for(i=0; i<nx; i++)
        for(j=0; j<ny; j++)
        {
          for(k=0; k<share; k++) 
            D->shareF[n][i][j][kstart+k]+=data[start+k]; 
          start+=share; 
        }
  }
  else if(rank%2==1 && rank!=D->N-1) {
    start=0; 
    for(n=0; n<numField; n++)
      for(i=0; i<nx; i++)
        for(j=0; j<ny; j++)
        {
          for(k=0; k<share; k++) {
            data[start+k]=D->shareF[n][i][j][kend+k];
            D->shareF[n][i][j][kend+k]=0.0;
          }
          start+=share; 
        }
    MPI_Send(data,num,MPI_DOUBLE,D->nextZrank,myrank,MPI_COMM_WORLD); 
  }            
  MPI_Barrier(MPI_COMM_WORLD);
  
  free(data);
}

void MPI_TransferJ_Period_Z(Domain *D,int numField)
{
  int n,i,j,k,start,numShare,nx,ny,share;
  int istart,iend,jstart,jend,kstart,kend;
  int myrank, nTasks,rank;
  double *data,tmp; 

  MPI_Status status;         

  istart=D->istart;    iend=D->iend;
  jstart=D->jstart;    jend=D->jend;
  kstart=D->kstart;    kend=D->kend;
  nx=1;	ny=1;
  if(D->dimension>1) { nx=iend+3; } else ;
  if(D->dimension>2) { ny=jend+3; } else ;
 
  MPI_Comm_size(MPI_COMM_WORLD, &nTasks);
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);            

  rank=(myrank%(D->M*D->N))/D->M;
  share=D->shareZ;

  // Bottom to Up 
	numShare=numField*nx*ny*(share-1);
  data = (double *)malloc(numShare*sizeof(double ));

  start=0; 
  for(n=0; n<numField; n++)
    for(i=0; i<nx; i++)
      for(j=0; j<ny; j++)
      {
        for(k=1; k<share; k++) {
          data[start+k-1]=D->shareF[n][i][j][kstart-k]; 
          D->shareF[n][i][j][kstart-k]=0.0;
        } 
        start+=share-1; 
      }

	if(D->N==1) {
    start=0; 
    for(n=0; n<numField; n++)
      for(i=0; i<nx; i++)
        for(j=0; j<ny; j++)
        {
          for(k=1; k<share; k++) 
            D->shareF[n][i][j][kend-k]+=data[start+k-1]; 
          start+=share-1; 
        }
	} else {
		if(rank==D->N-1)
    {
      MPI_Recv(data,numShare,MPI_DOUBLE,myrank-(D->N-1)*D->M,myrank-(D->N-1)*D->M, MPI_COMM_WORLD,&status);  
      start=0; 
      for(n=0; n<numField; n++)
        for(i=0; i<nx; i++)
          for(j=0; j<ny; j++)
          {
            for(k=1; k<share; k++) 
              D->shareF[n][i][j][kend-k]+=data[start+k-1]; 
            start+=share-1; 
          }
    } else if(rank==0)
      MPI_Send(data,numShare,MPI_DOUBLE,myrank+(D->N-1)*D->M,myrank,MPI_COMM_WORLD);             
    MPI_Barrier(MPI_COMM_WORLD);
	}
  free(data);

  // Up to Bottom
  numShare=numField*nx*ny*share;
  data = (double *)malloc(numShare*sizeof(double ));

  start=0; 
  for(n=0; n<numField; n++)
    for(i=0; i<nx; i++)
      for(j=0; j<ny; j++)
      {
        for(k=0; k<share; k++) {
          data[start+k]=D->shareF[n][i][j][kend+k]; 
          D->shareF[n][i][j][kend+k]=0.0;
        } 
        start+=share; 
      }
      
	if(D->N==1) {
    start=0;
    for(n=0; n<numField; n++)
      for(i=0; i<nx; i++)
        for(j=0; j<ny; j++)
        {
          for(k=0; k<share; k++) 
            D->shareF[n][i][j][kstart+k]+=data[start+k]; 
          start+=share; 
        }
	} else {
		if(rank==0)
    {
      MPI_Recv(data,numShare,MPI_DOUBLE,myrank+(D->N-1)*D->M,myrank+(D->N-1)*D->M, MPI_COMM_WORLD,&status);  
      start=0;
      for(n=0; n<numField; n++)
        for(i=0; i<nx; i++)
          for(j=0; j<ny; j++)
          {
            for(k=0; k<share; k++) 
              D->shareF[n][i][j][kstart+k]+=data[start+k]; 
            start+=share; 
          }
    } else if(rank==D->M-1)
      MPI_Send(data,numShare,MPI_DOUBLE,myrank-(D->N-1)*D->M,myrank, MPI_COMM_WORLD);             
    MPI_Barrier(MPI_COMM_WORLD);
	}
  free(data);
}
