#include <stdio.h>
#include <stdlib.h>
#include "mesh.h"
#include "mpi.h"

void MPI_TransferF_Yminus(Domain *D,int numField)
{
  int n,i,j,k,numberData,start,end,num,nx,nz,share;
  int istart,iend,jstart,jend,kstart,kend;
  int myrank, nTasks,rank; 
  double *data;

  MPI_Status status;         

  istart=D->istart;    iend=D->iend;
  jstart=D->jstart;    jend=D->jend;
  kstart=D->kstart;    kend=D->kend;

  MPI_Comm_size(MPI_COMM_WORLD, &nTasks);
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);            

  rank=(myrank%(D->M*D->N))%D->M;
  nx=1;	nz=1;
  if(D->dimension>1) { nx=iend+3; } else ;
  if(D->dimension>2) { nz=kend+3; } else ;

  share=D->shareY;
  num=numField*nx*nz*share;
  data = (double *)malloc(num*sizeof(double ));

  //Transferring even ~ odd cores 
	if(rank%2==0 && rank!=D->M-1)
  {
	  MPI_Recv(data,num,MPI_DOUBLE,D->nextYrank,D->nextYrank, MPI_COMM_WORLD,&status);  
	  start=0; 
	    //for(i=istart; i<iend; i++)
  	  for(n=0; n<numField; n++)	
	    for(i=0; i<nx; i++)
	      for(j=0; j<share; j++)	
		   {
		     for(k=0; k<nz; k++) 
             D->shareF[n][i][jend+j][k]=data[start+k]; 
           start+=nz; 
         } 
  }
  else if(rank%2==1)
  {
    start=0; 
	 for(n=0; n<numField; n++)	
	   for(i=0; i<nx; i++)
	     for(j=0; j<share; j++)	
	     {
          for(k=0; k<nz; k++) 
            data[start+k]=D->shareF[n][i][j+jstart][k]; 
          start+=nz; 
        }
    MPI_Send(data,num,MPI_DOUBLE,D->prevYrank,myrank,MPI_COMM_WORLD);             
  }
  MPI_Barrier(MPI_COMM_WORLD);

  //Transferring odd ~ even cores             
  if(rank%2==1 && rank!=D->M-1)
  {
	  MPI_Recv(data,num,MPI_DOUBLE,D->nextYrank,D->nextYrank, MPI_COMM_WORLD,&status);  
	  start=0; 
  	  for(n=0; n<numField; n++)	
	    for(i=0; i<nx; i++)
		   for(j=0; j<share; j++)	
			{
			  for(k=0; k<nz; k++) 
             D->shareF[n][i][jend+j][k]=data[start+k]; 
           start+=nz; 
         } 
  }
  else if(rank%2==0 && rank!=0)
  {
	 start=0; 
	 for(n=0; n<numField; n++)	
	   for(i=0; i<nx; i++)
	     for(j=0; j<share; j++)	
	     {
		    for(k=0; k<nz; k++) 
            data[start+k]=D->shareF[n][i][j+jstart][k]; 
          start+=nz; 
        }        
    MPI_Send(data,num,MPI_DOUBLE,D->prevYrank,myrank, MPI_COMM_WORLD);             
  }
  MPI_Barrier(MPI_COMM_WORLD);
  free(data);
}

void MPI_TransferF_Yplus(Domain *D,int numField)
{
  int n,i,j,k,num,nx,nz,share;
  int istart,iend,jstart,jend,kstart,kend;
  int myrank, nTasks,rank,start; 
  double *data;

  MPI_Status status;         
  MPI_Comm_size(MPI_COMM_WORLD, &nTasks);
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);            
  
  istart=D->istart;    iend=D->iend;
  jstart=D->jstart;    jend=D->jend;
  kstart=D->kstart;    kend=D->kend;

  rank=(myrank%(D->M*D->N))%D->M;
  nx=1;	nz=1;
  if(D->dimension>1) { nx=iend+3; } else ;
  if(D->dimension>2) { nz=kend+3; } else ;

  share=D->shareY;
	num=numField*nx*nz*(share-1);
  data = (double *)malloc(num*sizeof(double ));

  //Transferring even ~ odd cores 
  if(rank%2==1)
  {
	  MPI_Recv(data,num,MPI_DOUBLE,D->prevYrank,D->prevYrank, MPI_COMM_WORLD,&status);  
	  start=0;
     for(n=0; n<numField; n++)	
	    for(i=0; i<nx; i++)
		   for(j=1; j<share; j++)	
			{
			  for(k=0; k<nz; k++) 
             D->shareF[n][i][jstart-j][k]=data[start+k]; 
           start+=nz; 
         }
  }
  else if(rank%2==0 && rank!=D->M-1)
  {
    start=0; 
    for(n=0; n<numField; n++)	
	   for(i=0; i<nx; i++)
	     for(j=1; j<share; j++)	
	     {
		    for(k=0; k<nz; k++) 
            data[start+k]=D->shareF[n][i][jend-j][k]; 
          start+=nz; 
        }    
    MPI_Send(data,num,MPI_DOUBLE,D->nextYrank,myrank, MPI_COMM_WORLD);             
  }
  MPI_Barrier(MPI_COMM_WORLD);

  //Transferring odd ~ even cores             
  if(rank%2==0 && rank!=0)
  {
	  MPI_Recv(data,num,MPI_DOUBLE,D->prevYrank,D->prevYrank,MPI_COMM_WORLD,&status);  
	  start=0;
     for(n=0; n<numField; n++)	
	    for(i=0; i<nx; i++)
		   for(j=1; j<share; j++)	
			{
			  for(k=0; k<nz; k++) 
             D->shareF[n][i][jstart-j][k]=data[start+k]; 
           start+=nz; 
         }
  }
  else if(rank%2==1 && rank!=D->M-1)
  {
    start=0; 
    for(n=0; n<numField; n++)	
	   for(i=0; i<nx; i++)
	     for(j=1; j<share; j++)	
	     {
		    for(k=0; k<nz; k++) 
            data[start+k]=D->shareF[n][i][jend-j][k]; 
          start+=nz; 
        }      
    MPI_Send(data,num,MPI_DOUBLE,D->nextYrank,myrank,MPI_COMM_WORLD);             
  }
  MPI_Barrier(MPI_COMM_WORLD);
  free(data);
}

void MPI_TransferF_Period_Y(Domain *D,int numField)
{
  int n,i,j,k,numberData,start,end,num,nx,nz,share;
  int istart,iend,jstart,jend,kstart,kend;
  int myrank, nTasks,rank; 
  double *data;

  MPI_Status status;         

  istart=D->istart;    iend=D->iend;
  jstart=D->jstart;    jend=D->jend;
  kstart=D->kstart;    kend=D->kend;
  nx=1;	nz=1;
  if(D->dimension>1) { nx=iend+3; } else ;
  if(D->dimension>2) { nz=kend+3; } else ;    

  MPI_Comm_size(MPI_COMM_WORLD, &nTasks);
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);            

  rank=(myrank%(D->M*D->N))%D->M;
  share=D->shareY;

  // Bottom to Up 
  num=numField*nx*nz*share;
  data = (double *)malloc(num*sizeof(double ));

  start=0; 
  for(n=0; n<numField; n++)	
    for(i=0; i<nx; i++)
      for(j=0; j<share; j++)	
      {
        for(k=0; k<nz; k++) 
          data[start+k]=D->shareF[n][i][j+jstart][k]; 
        start+=nz; 
      }

	if(D->M==1) {
    start=0; 
    for(n=0; n<numField; n++)	
      for(i=0; i<nx; i++)
        for(j=0; j<share; j++)	
        {
          for(k=0; k<nz; k++) 
            D->shareF[n][i][jend+j][k]=data[start+k]; 
          start+=nz; 
        }
	} else {
		if(rank==0)
      MPI_Send(data,num,MPI_DOUBLE,myrank+(D->M-1),myrank,MPI_COMM_WORLD);             
    else if(rank==D->M-1) 	{
	    MPI_Recv(data,num,MPI_DOUBLE,myrank-(D->M-1),myrank-(D->M-1), MPI_COMM_WORLD,&status);  
   	  start=0; 
      for(n=0; n<numField; n++)	
        for(i=0; i<nx; i++)
          for(j=0; j<share; j++)	{
            for(k=0; k<nz; k++) 
              D->shareF[n][i][jend+j][k]=data[start+k]; 
            start+=nz; 
          }
    }
    MPI_Barrier(MPI_COMM_WORLD);
	}
	free(data);

  // Up to Bottom 
  num=numField*nx*nz*(share-1);
  data = (double *)malloc(num*sizeof(double ));

  start=0; 
  for(n=0; n<numField; n++)	
    for(i=0; i<nx; i++)
      for(j=1; j<share; j++)	
      {
        for(k=0; k<nz; k++) 
          data[start+k]=D->shareF[n][i][jend-j][k]; 
        start+=nz; 
      }

	if(D->M==1) {   
	  start=0;
     for(n=0; n<numField; n++)	
       for(i=0; i<nx; i++)
      	for(j=1; j<share; j++)	{
           for(k=0; k<nz; k++) 
             D->shareF[n][i][jstart-j][k]=data[start+k];
           start+=nz; 
         }
	} else {
		if(rank==D->M-1)
        MPI_Send(data,num,MPI_DOUBLE,myrank-(D->M-1),myrank, MPI_COMM_WORLD);             
      else if(rank==0)  {
      MPI_Recv(data,num,MPI_DOUBLE,myrank+(D->M-1),myrank+(D->M-1), MPI_COMM_WORLD,&status);  
      start=0;
      for(n=0; n<numField; n++)	
        for(i=0; i<nx; i++)
          for(j=1; j<share; j++)	
			 {
            for(k=0; k<nz; k++) 
              D->shareF[n][i][jstart-j][k]=data[start+k]; 
            start+=nz; 
          }
    }
    MPI_Barrier(MPI_COMM_WORLD);
	}
	free(data);
}



//---------------------------------------- Current or Density ---------------------------------------------//

void MPI_TransferJ_Yminus(Domain *D,int numField)
{
  int n,i,j,k,start,numShare,nx,nz,share;
  int istart,iend,jstart,jend,kstart,kend;
  int myrank, nTasks,rank;
  double *Yminus,tmp; 

  MPI_Status status;         

  istart=D->istart;    iend=D->iend;
  jstart=D->jstart;    jend=D->jend;
  kstart=D->kstart;    kend=D->kend;

  MPI_Comm_size(MPI_COMM_WORLD, &nTasks);
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);    
  nx=1;	nz=1;
  if(D->dimension>1) { nx=iend+3; } else ;
  if(D->dimension>2) { nz=kend+3; } else ;           

  share=D->shareY;
  rank=(myrank%(D->M*D->N))%D->M;
  numShare=numField*nx*nz*(share-1);
  Yminus = (double *)malloc(numShare*sizeof(double ));

  //Transferring even ~ odd cores 


  if(rank%2==0 && rank!=D->M-1)
  {
    MPI_Recv(Yminus,numShare,MPI_DOUBLE,D->nextYrank,D->nextYrank, MPI_COMM_WORLD,&status);  
    start=0; 
    for(n=0; n<numField; n++)
      for(i=0; i<nx; i++)
        for(j=1; j<share; j++)	
        {
          for(k=0; k<nz; k++) 
		      D->shareF[n][i][jend-j][k]+=Yminus[start+k]; 
		    start+=nz; 
		  }
  }
  else if(rank%2==1) 
  {
    start=0; 
    for(n=0; n<numField; n++)
      for(i=0; i<nx; i++)
        for(j=1; j<share; j++)	
        {
          for(k=0; k<nz; k++) {
	  	      Yminus[start+k]=D->shareF[n][i][jstart-j][k]; 
            D->shareF[n][i][jstart-j][k]=0.0;
          } 
	       start+=nz; 
	  	  }
    MPI_Send(Yminus,numShare,MPI_DOUBLE,D->prevYrank,myrank,MPI_COMM_WORLD);    
  }         
  MPI_Barrier(MPI_COMM_WORLD);

  //Transferring odd ~ even cores                     
  if(rank%2==1 && rank!=D->M-1)
  {
    MPI_Recv(Yminus,numShare,MPI_DOUBLE,D->nextYrank,D->nextYrank, MPI_COMM_WORLD,&status);  
    start=0; 
    for(n=0; n<numField; n++)
      for(i=0; i<nx; i++)
        for(j=1; j<share; j++)	
        {
          for(k=0; k<nz; k++) 
		      D->shareF[n][i][jend-j][k]+=Yminus[start+k]; 
			 start+=nz; 
		  }
  }
  else if(rank%2==0 && rank!=0) 
  {
    start=0; 
    for(n=0; n<numField; n++)
      for(i=0; i<nx; i++)
        for(j=1; j<share; j++)	
        {
          for(k=0; k<nz; k++) {
	  	      Yminus[start+k]=D->shareF[n][i][jstart-j][k]; 
            D->shareF[n][i][jstart-j][k]=0.0;
          } 
	       start+=nz; 
	  	  }
    MPI_Send(Yminus,numShare,MPI_DOUBLE,D->prevYrank,myrank, MPI_COMM_WORLD);    
  }         
  MPI_Barrier(MPI_COMM_WORLD);

  free(Yminus);
}

void MPI_TransferJ_Yplus(Domain *D,int numField)
{
  int n,i,j,k,start,numShare,nx,nz,share;
  int istart,iend,jstart,jend,kstart,kend;
  int myrank, nTasks,rank; 
  double *Yplus,tmp;

  MPI_Status status;         
   
  istart=D->istart;    iend=D->iend;
  jstart=D->jstart;    jend=D->jend;
  kstart=D->kstart;    kend=D->kend;
  nx=1;	nz=1;
  if(D->dimension>1) { nx=iend+3; } else ;
  if(D->dimension>2) { nz=kend+3; } else ;    

  MPI_Comm_size(MPI_COMM_WORLD, &nTasks);
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);            

  share=D->shareY;
  rank=(myrank%(D->M*D->N))%D->M;
  numShare=numField*nx*nz*share;
  Yplus = (double *)malloc(numShare*sizeof(double ));

  //Transferring even ~ odd cores 

      
  if(rank%2==1)
  {
    MPI_Recv(Yplus,numShare,MPI_DOUBLE,D->prevYrank,D->prevYrank, MPI_COMM_WORLD,&status);  
    start=0;
    for(n=0; n<numField; n++)
      for(i=0; i<nx; i++)
        for(j=0; j<share; j++)	
		  {
          for(k=0; k<nz; k++) 
	         D->shareF[n][i][jstart+j][k]+=Yplus[start+k]; 
		    start+=nz; 
	     } 
  }
  else if(rank%2==0 && rank!=D->M-1) 
  {
    start=0; 
    for(n=0; n<numField; n++)
      for(i=0; i<nx; i++)
        for(j=0; j<share; j++)	
		  {
          for(k=0; k<nz; k++) {
	  	      Yplus[start+k]=D->shareF[n][i][jend+j][k]; 
            D->shareF[n][i][jend+j][k]=0.0;
          } 
	  	    start+=nz; 
	  	  }
    MPI_Send(Yplus,numShare,MPI_DOUBLE,D->nextYrank,myrank, MPI_COMM_WORLD);  
  }           
  MPI_Barrier(MPI_COMM_WORLD);

  //Transferring odd ~ even cores             
  if(rank%2==0 && rank!=0)
  {
    MPI_Recv(Yplus,numShare,MPI_DOUBLE,D->prevYrank,D->prevYrank,MPI_COMM_WORLD,&status);  
    start=0;
    for(n=0; n<numField; n++)
      for(i=0; i<nx; i++)
        for(j=0; j<share; j++)	
        {
          for(k=0; k<nz; k++) 
	         D->shareF[n][i][jstart+j][k]+=Yplus[start+k]; 
		    start+=nz; 
	     } 
  }
  else if(rank%2==1 && rank!=D->M-1) {
    start=0; 
    for(n=0; n<numField; n++)
      for(i=0; i<nx; i++)
        for(j=0; j<share; j++)	
        {
          for(k=0; k<nz; k++) {
	  	      Yplus[start+k]=D->shareF[n][i][jend+j][k]; 
            D->shareF[n][i][jend+j][k]=0.0;
          } 
	  	    start+=nz; 
	  	  }
    MPI_Send(Yplus,numShare,MPI_DOUBLE,D->nextYrank,myrank,MPI_COMM_WORLD);   
  }          
  MPI_Barrier(MPI_COMM_WORLD);

  free(Yplus);
}


void MPI_TransferJ_Period_Y(Domain *D,int numField)
{
  int n,i,j,k,start,numShare,nx,nz,share;
  int istart,iend,jstart,jend,kstart,kend;
  int myrank, nTasks,rank;
  double *data,tmp; 

  MPI_Status status;         

  istart=D->istart;    iend=D->iend;
  jstart=D->jstart;    jend=D->jend;
  kstart=D->kstart;    kend=D->kend;
  nx=1;	nz=1;
  if(D->dimension>1) { nx=iend+3; } else ;
  if(D->dimension>2) { nz=kend+3; } else ;    
 
  MPI_Comm_size(MPI_COMM_WORLD, &nTasks);
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);            

  rank=(myrank%(D->M*D->N))%D->M;
  share=D->shareY;
  
  // Bottom to Up 
	numShare=numField*nx*nz*(share-1);
  data = (double *)malloc(numShare*sizeof(double ));

  start=0; 
  for(n=0; n<numField; n++)
    for(i=0; i<nx; i++)
      for(j=1; j<share; j++)	
      {
        for(k=0; k<nz; k++) {
          data[start+k]=D->shareF[n][i][jstart-j][k]; 
          D->shareF[n][i][jstart-j][k]=0.0;
        } 
        start+=nz; 
      }

	if(D->M==1) {
    start=0; 
    for(n=0; n<numField; n++)
      for(i=0; i<nx; i++)
        for(j=1; j<share; j++)	
        { 
          for(k=0; k<nz; k++) 
            D->shareF[n][i][jend-j][k]+=data[start+k]; 
          start+=nz; 
        }
	} else {
		if(rank==D->M-1)
    {
      MPI_Recv(data,numShare,MPI_DOUBLE,myrank-(D->M-1),myrank-(D->M-1), MPI_COMM_WORLD,&status);  
      start=0; 
      for(n=0; n<numField; n++)
        for(i=0; i<nx; i++)
          for(j=1; j<share; j++)	
			 {
            for(k=0; k<nz; k++) 
              D->shareF[n][i][jend-j][k]+=data[start+k]; 
            start+=nz; 
          }
    } else if(rank==0)
      MPI_Send(data,numShare,MPI_DOUBLE,myrank+(D->M-1),myrank,MPI_COMM_WORLD);             
    MPI_Barrier(MPI_COMM_WORLD);
	}
  free(data);

  // Up to Bottom
  numShare=numField*nx*nz*share;
  data = (double *)malloc(numShare*sizeof(double ));

  start=0; 
  for(n=0; n<numField; n++)
    for(i=0; i<nx; i++)
      for(j=0; j<share; j++)	
      {
        for(k=0; k<nz; k++) {
          data[start+k]=D->shareF[n][i][jend+j][k]; 
          D->shareF[n][i][jend+j][k]=0.0;
        } 
        start+=nz; 
      }
      
	if(D->M==1) {
     start=0;
     for(n=0; n<numField; n++)
       for(i=0; i<nx; i++)
         for(j=0; j<share; j++)	
         {
           for(k=0; k<nz; k++) 
             D->shareF[n][i][jstart+j][k]+=data[start+k]; 
           start+=nz; 
         }
	} else 	{
	  if(rank==0)
     {
       MPI_Recv(data,numShare,MPI_DOUBLE,myrank+(D->M-1),myrank+(D->M-1), MPI_COMM_WORLD,&status);  
       start=0;
       for(n=0; n<numField; n++)
         for(i=0; i<nx; i++)
           for(j=0; j<share; j++)	
			  {
             for(k=0; k<nz; k++) 
               D->shareF[n][i][jstart+j][k]+=data[start+k]; 
             start+=nz; 
           }
     } else if(rank==D->M-1)
       MPI_Send(data,numShare,MPI_DOUBLE,myrank-(D->M-1),myrank, MPI_COMM_WORLD);             
     MPI_Barrier(MPI_COMM_WORLD);
	}
   free(data);
}

