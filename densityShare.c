#include <stdio.h>
#include <stdlib.h>
#include "mesh.h"
#include "mpi.h"

void MPI_Density_Zminus(Domain *D,int nx,int ny,int share)
{
  int i,j,k,numberData,start,end,ibegin,ibottom;
  int istart,iend,jstart,jend,kstart,kend;
  int myrank, nTasks,rank; 

  MPI_Status status;         
  MPI_Comm_size(MPI_COMM_WORLD, &nTasks);
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);            

  ibegin=0;
  ibottom=nx;   

  istart=D->istart;
  iend=D->iend;
  jstart=D->jstart;
  jend=D->jend;
  kstart=D->kstart;
  kend=D->kend;
  //first 3 is 3 layer, 2nd 3 is 6 field variables
//    numberData=share*6*nx*nz; 
  rank=(int)(myrank/D->M);   

  //Transferring even ~ odd cores 
  start=0; 
  for(k=1; k<share; k++)	// 3 layers in y direction
    for(j=0; j<ny; j++)
    {
      for(i=ibegin; i<ibottom; i++)
        D->minusDenZ[start+i]=D->Rho[i][j][kstart-k];
      start+=nx;
    }

  if(rank%2==0 && rank!=D->N-1)
  {
    MPI_Recv(D->minusDenZ,D->numMinusDenZ, MPI_DOUBLE, myrank+D->M, myrank+D->M, MPI_COMM_WORLD,&status);  
    start=0; 
    for(k=1; k<share; k++)	// 3 layers in y direction
      for(j=0; j<ny; j++)
      {
        for(i=ibegin; i<ibottom; i++)
          D->Rho[i][j][kend-1-k]+=D->minusDenZ[start+i];
        start+=nx;
      }  
  }
  else if(rank%2==1)
    MPI_Send(D->minusDenZ,D->numMinusDenZ, MPI_DOUBLE, myrank-D->M, myrank, MPI_COMM_WORLD);             
  MPI_Barrier(MPI_COMM_WORLD);

  //Transferring odd ~ even cores             
  start=0; 
  for(k=1; k<share; k++)	// 3 layers in y direction
    for(j=0; j<ny; j++)
    {
      for(i=ibegin; i<ibottom; i++)
        D->minusDenZ[start+i]=D->Rho[i][j][kstart-k];
      start+=nx;
    }
        
  if(rank%2==1 && rank!=D->N-1)
  {
    MPI_Recv(D->minusDenZ,D->numMinusDenZ, MPI_DOUBLE, myrank+D->M, myrank+D->M, MPI_COMM_WORLD,&status);  
    start=0; 
    for(k=1; k<share; k++)	// 1 layers in y direction
      for(j=0; j<ny; j++)
      {
        for(i=ibegin; i<ibottom; i++)
          D->Rho[i][j][kend-1-k]+=D->minusDenZ[start+i];
        start+=nx;
      }  
  }
  else if(rank%2==0 && rank!=0)
    MPI_Send(D->minusDenZ,D->numMinusDenZ, MPI_DOUBLE, myrank-D->M, myrank, MPI_COMM_WORLD);             
  MPI_Barrier(MPI_COMM_WORLD);
}


void MPI_Density_Zplus(Domain *D,int nx,int ny,int share)
{
  int i,j,k,start,end,ibegin,ibottom;
  int istart,iend,jstart,jend,kstart,kend;
  int myrank, nTasks,rank; 

  MPI_Status status;         
  MPI_Comm_size(MPI_COMM_WORLD, &nTasks);
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);            
   
  ibegin=0;
  ibottom=nx;

  istart=D->istart;
  iend=D->iend;
  jstart=D->jstart;
  jend=D->jend;
  kstart=D->kstart;
  kend=D->kend;
  //first 2 is 2 layer, 2nd 6 is 6 field variables
//    numberData=(share-1)*6*nx*nz; 
  rank=(int)(myrank/D->M);   

  //Transferring even ~ odd cores 
  start=0; 
  for(k=0; k<share; k++)
    for(j=0; j<ny; j++)
    {
      for(i=ibegin; i<ibottom; i++)
        D->plusDenZ[start+i]=D->Rho[i][j][kend+k];
      start+=nx;
    }
      
  if(rank%2==1)
  {
     MPI_Recv(D->plusDenZ,D->numPlusDenZ, MPI_DOUBLE, myrank-D->M, myrank-D->M, MPI_COMM_WORLD,&status);  
     start=0;
     for(k=0; k<share; k++)
       for(j=0; j<ny; j++)
       {
         for(i=ibegin; i<ibottom; i++)
           D->Rho[i][j][kstart-k]+=D->plusDenZ[start+i];
         start+=nx;
       }
  }
  else if(rank%2==0 && rank!=D->N-1)
    MPI_Send(D->plusDenZ,D->numPlusDenZ, MPI_DOUBLE, myrank+D->M, myrank, MPI_COMM_WORLD);             
  MPI_Barrier(MPI_COMM_WORLD);

  //Transferring odd ~ even cores             
  start=0; 
  for(k=0; k<share; k++)
    for(j=0; j<ny; j++)
    {
      for(i=ibegin; i<ibottom; i++)
        D->plusDenZ[start+i]=D->Rho[i][j][kend+k];
      start+=nx;
    }
        
  if(rank%2==0 && rank!=0)
  {
    MPI_Recv(D->plusDenZ,D->numPlusDenZ, MPI_DOUBLE, myrank-D->M, myrank-D->M, MPI_COMM_WORLD,&status);  
    start=0;
    for(k=0; k<share; k++)
      for(j=0; j<ny; j++)
      {
        for(i=ibegin; i<ibottom; i++)
          D->Rho[i][j][kstart+k]+=D->plusDenZ[start+i];
        start+=nx;
      }
  }
  else if(rank%2==1 && rank!=D->N-1)
    MPI_Send(D->plusDenZ,D->numPlusDenZ, MPI_DOUBLE, myrank+D->M, myrank, MPI_COMM_WORLD);             
  MPI_Barrier(MPI_COMM_WORLD);
}

void MPI_Density_Yminus(Domain *D,int nx,int nz,int share)
{
  int i,j,k,numberData,start,end,ibegin,ibottom;
  int istart,iend,jstart,jend,kstart,kend;
  int myrank, nTasks,rank; 

  MPI_Status status;         
  MPI_Comm_size(MPI_COMM_WORLD, &nTasks);
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);            

  ibegin=0;
  ibottom=nx;   

  istart=D->istart;
  iend=D->iend;
  jstart=D->jstart;
  jend=D->jend;
  kstart=D->kstart;
  kend=D->kend;
  //first 3 is 3 layer, 2nd 3 is 6 field variables
//    numberData=share*6*nx*nz; 
  rank=myrank%D->M;   

  //Transferring even ~ odd cores 
  start=0; 
  for(j=1; j<share; j++)	// 3 layers in y direction
    for(k=0; k<nz; k++)
    {
      for(i=ibegin; i<ibottom; i++)
        D->minusDenY[start+i]=D->Rho[i][jstart-j][k];
      start+=nx;
    }

  if(rank%2==0 && rank!=D->M-1)
  {
    MPI_Recv(D->minusDenY,D->numMinusDenY, MPI_DOUBLE, myrank+1, myrank+1, MPI_COMM_WORLD,&status);  
    start=0; 
    for(j=1; j<share; j++)	// 3 layers in y direction
      for(k=0; k<nz; k++)
      {
        for(i=ibegin; i<ibottom; i++)
          D->Rho[i][jend-1-j][k]+=D->minusDenY[start+i];
        start+=nx;
      }  
  }
  else if(rank%2==1)
    MPI_Send(D->minusDenY,D->numMinusDenY, MPI_DOUBLE, myrank-1, myrank, MPI_COMM_WORLD);             
  MPI_Barrier(MPI_COMM_WORLD);

  //Transferring odd ~ even cores             
  start=0; 
  for(j=1; j<share; j++)	// 3 layers in y direction
    for(k=0; k<nz; k++)
    {
      for(i=ibegin; i<ibottom; i++)
        D->minusDenY[start+i]=D->Rho[i][jstart-j][k];
      start+=nx;
    }
        
  if(rank%2==1 && rank!=D->M-1)
  {
    MPI_Recv(D->minusDenY,D->numMinusDenY, MPI_DOUBLE, myrank+1, myrank+1, MPI_COMM_WORLD,&status);  
    start=0; 
    for(j=1; j<share; j++)	// 1 layers in y direction
      for(k=0; k<nz; k++)
      {
        for(i=ibegin; i<ibottom; i++)
          D->Rho[i][jend-1-j][k]+=D->minusDenY[start+i];
        start+=nx;
      }  
  }
  else if(rank%2==0 && rank!=0)
    MPI_Send(D->minusDenY,D->numMinusDenY, MPI_DOUBLE, myrank-1, myrank, MPI_COMM_WORLD);             
  MPI_Barrier(MPI_COMM_WORLD);
}


void MPI_Density_Yplus(Domain *D,int nx,int nz,int share)
{
  int i,j,k,start,end,ibegin,ibottom;
  int istart,iend,jstart,jend,kstart,kend;
  int myrank, nTasks,rank; 

  MPI_Status status;         
  MPI_Comm_size(MPI_COMM_WORLD, &nTasks);
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);            
   
  ibegin=0;
  ibottom=nx;

  istart=D->istart;
  iend=D->iend;
  jstart=D->jstart;
  jend=D->jend;
  kstart=D->kstart;
  kend=D->kend;
  //first 2 is 2 layer, 2nd 6 is 6 field variables
//    numberData=(share-1)*6*nx*nz; 
  rank=myrank%D->M;   

  //Transferring even ~ odd cores 
  start=0; 
  for(j=0; j<share; j++)
    for(k=0; k<nz; k++)
    {
      for(i=ibegin; i<ibottom; i++)
        D->plusDenY[start+i]=D->Rho[i][jend+j][k];
      start+=nx;
    }
      
  if(rank%2==1)
  {
     MPI_Recv(D->plusDenY,D->numPlusDenY, MPI_DOUBLE, myrank-1, myrank-1, MPI_COMM_WORLD,&status);  
     start=0;
     for(j=0; j<share; j++)
       for(k=0; k<nz; k++)
       {
         for(i=ibegin; i<ibottom; i++)
           D->Rho[i][jstart-j][k]+=D->plusDenY[start+i];
         start+=nx;
       }
  }
  else if(rank%2==0 && rank!=D->M-1)
    MPI_Send(D->plusDenY,D->numPlusDenY, MPI_DOUBLE, myrank+1, myrank, MPI_COMM_WORLD);             
  MPI_Barrier(MPI_COMM_WORLD);

  //Transferring odd ~ even cores             
  start=0; 
  for(j=0; j<share; j++)
    for(k=0; k<nz; k++)
    {
      for(i=ibegin; i<ibottom; i++)
        D->plusDenY[start+i]=D->Rho[i][jend+j][k];
      start+=nx;
    }
        
  if(rank%2==0 && rank!=0)
  {
    MPI_Recv(D->plusDenY,D->numPlusDenY, MPI_DOUBLE, myrank-1, myrank-1, MPI_COMM_WORLD,&status);  
    start=0;
    for(j=0; j<share; j++)
      for(k=0; k<nz; k++)
      {
        for(i=ibegin; i<ibottom; i++)
          D->Rho[i][jstart+j][k]+=D->plusDenY[start+i];
        start+=nx;
      }
  }
  else if(rank%2==1 && rank!=D->M-1)
    MPI_Send(D->plusDenY,D->numPlusDenY, MPI_DOUBLE, myrank+1, myrank, MPI_COMM_WORLD);             
  MPI_Barrier(MPI_COMM_WORLD);
}
