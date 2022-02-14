#include <stdio.h>
#include <stdlib.h>
#include "mesh.h"
#include "mpi.h"


void MPI_TransferF_Xminus(Domain *D,int numField)
{
  int n,i,j,k,num,start,end,ny,nz,share;
  int istart,iend,jstart,jend,kstart,kend;
  int myrank, nTasks,rank; 
  double *data;

  MPI_Status status;         

  istart=D->istart;    iend=D->iend;
  jstart=D->jstart;    jend=D->jend;
  kstart=D->kstart;    kend=D->kend;

  MPI_Comm_size(MPI_COMM_WORLD, &nTasks);
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);            

  rank=myrank/(D->M*D->N);
	ny=1;	nz=1;
	if(D->dimension>1) {	ny=jend+3; } else ;
	if(D->dimension>2) {	nz=kend+3; } else ;
  share=D->shareX;

  num=numField*ny*nz*share;
  data = (double *)malloc(num*sizeof(double ));
 
  //Transferring even ~ odd cores 
  start=0; 
  for(n=0; n<numField; n++)
    for(i=0; i<share; i++)
      for(j=0; j<ny; j++)   {
	      for(k=0; k<nz; k++) 
          data[start+k]=D->shareF[n][i+istart][j][k]; 
        start+=nz; 
      }

  if(rank%2==0 && rank!=D->L-1)
  {
    MPI_Recv(data,num,MPI_DOUBLE,D->nextXrank,D->nextXrank, MPI_COMM_WORLD,&status);  
    start=0; 
    for(n=0; n<numField; n++)
      for(i=0; i<share; i++)
        for(j=0; j<ny; j++)  {
          for(k=0; k<nz; k++) 
            D->shareF[n][iend+i][j][k]=data[start+k]; 
          start+=nz; 
        }
      }
  else if(rank%2==1)
    MPI_Send(data,num,MPI_DOUBLE,D->prevXrank,myrank,MPI_COMM_WORLD);             
  MPI_Barrier(MPI_COMM_WORLD);

  //Transferring odd ~ even cores             
  start=0; 
  for(n=0; n<numField; n++)
    for(i=0; i<share; i++)
      for(j=0; j<ny; j++)   {
	      for(k=0; k<nz; k++) 
          data[start+k]=D->shareF[n][i+istart][j][k]; 
        start+=nz; 
      }

        
  if(rank%2==1 && rank!=D->L-1)
  {
    MPI_Recv(data,num,MPI_DOUBLE,D->nextXrank,D->nextXrank, MPI_COMM_WORLD,&status);  
    start=0; 
    for(n=0; n<numField; n++)
      for(i=0; i<share; i++)
        for(j=0; j<ny; j++)  {
          for(k=0; k<nz; k++) 
            D->shareF[n][iend+i][j][k]=data[start+k]; 
          start+=nz; 
        }
  }
  else if(rank%2==0 && rank!=0)
    MPI_Send(data,num,MPI_DOUBLE,D->prevXrank,myrank, MPI_COMM_WORLD);             
  MPI_Barrier(MPI_COMM_WORLD);
  free(data);
}

void MPI_TransferF_Xplus(Domain *D,int numField)
{
  int n,i,j,k,num,ny,nz,share;
  int istart,iend,jstart,jend,kstart,kend;
  int myrank, nTasks,rank,start; 
  double *data;

  MPI_Status status;         
  
  istart=D->istart;   iend=D->iend;
  jstart=D->jstart;   jend=D->jend;
  kstart=D->kstart;   kend=D->kend;

  MPI_Comm_size(MPI_COMM_WORLD, &nTasks);
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);            

  rank=myrank/(D->M*D->N);
	ny=1;	nz=1;
	if(D->dimension>1) {	ny=jend+3; } else ;
	if(D->dimension>2) {	nz=kend+3; } else ;
  share=D->shareX;
  num=numField*ny*nz*2;
  data = (double *)malloc(num*sizeof(double ));

  //Transferring even ~ odd cores 
  start=0; 
  for(n=0; n<numField; n++)
    for(i=1; i<share; i++)
      for(j=0; j<ny; j++)	{
        for(k=0; k<nz; k++) 
          data[start+k]=D->shareF[n][iend-i][j][k]; 
        start+=nz;
      }       
     
  if(rank%2==1)
  {
    MPI_Recv(data,num,MPI_DOUBLE,D->prevXrank,D->prevXrank, MPI_COMM_WORLD,&status);  
    start=0;
    for(n=0; n<numField; n++)
      for(i=1; i<share; i++)
        for(j=0; j<ny; j++)	{
          for(k=0; k<nz; k++) 
            D->shareF[n][istart-i][j][k]=data[start+k]; 
          start+=nz; 
        }
    }
  else if(rank%2==0 && rank!=D->L-1)
     MPI_Send(data,num,MPI_DOUBLE,D->nextXrank,myrank, MPI_COMM_WORLD);             
	MPI_Barrier(MPI_COMM_WORLD);

  start=0; 
  for(n=0; n<numField; n++)
    for(i=1; i<share; i++)
      for(j=0; j<ny; j++)	{
        for(k=0; k<nz; k++) 
          data[start+k]=D->shareF[n][iend-i][j][k]; 
        start+=nz;
      }      
									 
  if(rank%2==0 && rank!=0)
  {
	  MPI_Recv(data,num,MPI_DOUBLE,D->prevXrank,D->prevXrank,MPI_COMM_WORLD,&status);  
    start=0;
    for(n=0; n<numField; n++)
      for(i=1; i<share; i++)
        for(j=0; j<ny; j++)	{
          for(k=0; k<nz; k++) 
            D->shareF[n][istart-i][j][k]=data[start+k]; 
          start+=nz; 
        }
  }
  else if(rank%2==1 && rank!=D->L-1)
    MPI_Send(data,num,MPI_DOUBLE,D->nextXrank,myrank,MPI_COMM_WORLD);             
   
  free(data);
}

void MPI_TransferF_Period_X(Domain *D,int numField)
{
  int n,i,j,k,num,start,end,ny,nz,share;
  int istart,iend,jstart,jend,kstart,kend;
  int myrank, nTasks,rank; 
  double *data;

  MPI_Status status;         

  istart=D->istart;    iend=D->iend;
  jstart=D->jstart;    jend=D->jend;
  kstart=D->kstart;    kend=D->kend;
	ny=1;	nz=1;
	if(D->dimension>1) {	ny=jend+3; } else ;
	if(D->dimension>2) {	nz=kend+3; } else ;

  MPI_Comm_size(MPI_COMM_WORLD, &nTasks);
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);            
  share=D->shareX;
  rank=myrank/(D->M*D->N);

  // Left to Right 
  num=numField*ny*nz*share;
  data = (double *)malloc(num*sizeof(double ));
 
  start=0; 
  for(n=0; n<numField; n++)
    for(i=0; i<share; i++)
      for(j=0; j<ny; j++)   {
        for(k=0; k<nz; k++) 
          data[start+k]=D->shareF[n][i+istart][j][k]; 
        start+=nz; 
      }

	if(D->L==1) {
		start=0; 
    for(n=0; n<numField; n++)
      for(i=0; i<share; i++)
        for(j=0; j<ny; j++)        {
          for(k=0; k<nz; k++) 
            D->shareF[n][iend+i][j][k]=data[start+k]; 
          start+=nz; 
        }
	} else {
		if(rank==D->L-1)   {
      MPI_Recv(data,num,MPI_DOUBLE,myrank-(D->L-1)*D->M*D->N,myrank-(D->L-1)*D->M*D->N, MPI_COMM_WORLD,&status);  
      start=0; 
      for(n=0; n<numField; n++)
        for(i=0; i<share; i++)
          for(j=0; j<ny; j++)        {
            for(k=0; k<nz; k++) 
              D->shareF[n][iend+i][j][k]=data[start+k]; 
            start+=nz; 
          }
    }  else if(rank==0)
      MPI_Send(data,num,MPI_DOUBLE,myrank+(D->L-1)*D->M*D->N,myrank,MPI_COMM_WORLD);             
    MPI_Barrier(MPI_COMM_WORLD);
	}
	free(data);

  // Rgitht to Left 
  num=numField*ny*nz*(share-1);
  data = (double *)malloc(num*sizeof(double ));

  start=0; 
  for(n=0; n<numField; n++)
    for(i=1; i<share; i++)
      for(j=0; j<ny; j++)	{
        for(k=0; k<nz; k++) 
          data[start+k]=D->shareF[n][iend-i][j][k]; 
        start+=nz; 
      }

  if(D->L==1) {
		start=0;
    for(n=0; n<numField; n++)
      for(i=1; i<share; i++)
      	for(j=0; j<ny; j++)	{
          for(k=0; k<nz; k++) 
            D->shareF[n][istart-i][j][k]=data[start+k]; 
          start+=nz; 
        }
	} else {      
		if(rank==0)   {
			MPI_Recv(data,num,MPI_DOUBLE,myrank+(D->L-1)*D->M*D->N,myrank+(D->L-1)*D->M*D->N, MPI_COMM_WORLD,&status);  
      start=0;
      for(n=0; n<numField; n++)
        for(i=1; i<share; i++)
          for(j=0; j<ny; j++)	{
            for(k=0; k<nz; k++) 
              D->shareF[n][istart-i][j][k]=data[start+k]; 
            start+=nz; 
          }
    }  else if(rank==D->L-1)
      MPI_Send(data,num,MPI_DOUBLE,myrank-(D->L-1)*D->M*D->N,myrank, MPI_COMM_WORLD);             
   	MPI_Barrier(MPI_COMM_WORLD);
	}
  
  free(data);
}



// ----------------------------------- Current, Density ----------------------------------------//

void MPI_TransferJ_Xplus(Domain *D,int numField)
{
  int n,i,j,k,start,ny,nz,share;
  int istart,iend,jstart,jend,kstart,kend;
  int myrank, nTasks, rank,num;
  double *data;

  MPI_Status status;

  MPI_Comm_size(MPI_COMM_WORLD, &nTasks);
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

  istart=D->istart;    iend=D->iend;
  jstart=D->jstart;    jend=D->jend;
  kstart=D->kstart;    kend=D->kend;
	ny=1;	nz=1;
	if(D->dimension>1) {	ny=jend+3; } else ;
	if(D->dimension>2) {	nz=kend+3; } else ;
  share=D->shareX;

  rank=myrank/(D->M*D->N);
  num=numField*ny*nz*share;
  data = (double *)malloc(num*sizeof(double ));

  //Transferring even ~ odd cores 


  if(rank%2==1)
  {
    MPI_Recv(data,num,MPI_DOUBLE,D->prevXrank,D->prevXrank,MPI_COMM_WORLD,&status);
    start=0;
    for(n=0; n<numField; n++)
      for(i=0; i<share; i++)
        for(j=0; j<ny; j++) {
          for(k=0; k<nz; k++) 
		        D->shareF[n][istart+i][j][k]+=data[k+start]; 
			    start+=nz; 
		    } 
  }
  else if(rank%2==0 && rank!=D->L-1) {
    start=0;
    for(n=0; n<numField; n++)
      for(i=0; i<share; i++)
        for(j=0; j<ny; j++)   
		  {
          for(k=0; k<nz; k++) {
		    	data[start+k]=D->shareF[n][iend+i][j][k]; 
            D->shareF[n][iend+i][j][k]=0.0;
          } 
		    start+=nz; 
		  }
    MPI_Send(data,num,MPI_DOUBLE,D->nextXrank,myrank, MPI_COMM_WORLD);
  }
  MPI_Barrier(MPI_COMM_WORLD);



  if(rank%2==0 && rank!=0)
  {
    MPI_Recv(data,num,MPI_DOUBLE,D->prevXrank,D->prevXrank,MPI_COMM_WORLD,&status);
    start=0;
    for(n=0; n<numField; n++)
      for(i=0; i<share; i++)
        for(j=0; j<ny; j++) {
          for(k=0; k<nz; k++) 
		        D->shareF[n][istart+i][j][k]+=data[k+start]; 
			    start+=nz; 
		    } 
  }
  else if(rank%2==1 && rank!=D->L-1) 
  {
    start=0;
    for(n=0; n<numField; n++)
      for(i=0; i<share; i++)
        for(j=0; j<ny; j++)  {
          for(k=0; k<nz; k++) {
	  	      data[start+k]=D->shareF[n][iend+i][j][k]; 
            D->shareF[n][iend+i][j][k]=0.0;
          } 
	  	    start+=nz; 
	  	  }
    MPI_Send(data,num,MPI_DOUBLE,D->nextXrank,myrank, MPI_COMM_WORLD);
  }
  MPI_Barrier(MPI_COMM_WORLD);
  free(data);
}
                                                

void MPI_TransferJ_Xminus(Domain *D,int numField)
{
  int n,i,j,k,start,ny,nz,share;
  int istart,iend,jstart,jend,kstart,kend;
  int myrank, nTasks, rank,num;
  double *data;

  MPI_Status status;

  MPI_Comm_size(MPI_COMM_WORLD, &nTasks);
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

  istart=D->istart;    iend=D->iend;
  jstart=D->jstart;    jend=D->jend;
  kstart=D->kstart;    kend=D->kend;
	ny=1;	nz=1;
	if(D->dimension>1) {	ny=jend+3; } else ;
	if(D->dimension>2) {	nz=kend+3; } else ;
  share=D->shareX;

  rank=myrank/(D->M*D->N);
  num=numField*ny*nz*(share-1);
  data = (double *)malloc(num*sizeof(double ));

  //Transferring even ~ odd cores 


  if(rank%2==0 && rank!=D->L-1)
  {
    MPI_Recv(data,num,MPI_DOUBLE,D->nextXrank,D->nextXrank,MPI_COMM_WORLD,&status);
    start=0;
    for(n=0; n<numField; n++)
      for(i=1; i<share; i++)
        for(j=0; j<ny; j++) {
          for(k=0; k<nz; k++) 
		        D->shareF[n][iend-i][j][k]+=data[k+start]; 
			    start+=nz; 
		    } 
  }
  else if(rank%2==1) 
  {
    start=0;
    for(n=0; n<numField; n++)
      for(i=1; i<share; i++)
        for(j=0; j<ny; j++) 
        {
          for(k=0; k<nz; k++) {
	  	      data[start+k]=D->shareF[n][istart-i][j][k]; 
            D->shareF[n][istart-i][j][k]=0.0;
          } 
	  	    start+=nz; 
	  	  }    
    MPI_Send(data,num,MPI_DOUBLE,D->prevXrank,myrank, MPI_COMM_WORLD);
  }
  MPI_Barrier(MPI_COMM_WORLD);

    //Transferring odd ~ even cores             


  if(rank%2==1 && rank!=D->L-1)
  {
    MPI_Recv(data,num,MPI_DOUBLE,D->nextXrank,D->nextXrank,MPI_COMM_WORLD,&status);
    start=0;
    for(n=0; n<numField; n++)
      for(i=1; i<share; i++)
        for(j=0; j<ny; j++) {
          for(k=0; k<nz; k++) 
		        D->shareF[n][iend-i][j][k]+=data[k+start]; 
			    start+=nz; 
		    } 
  }
  else if(rank%2==0 && rank!=0) 
  {
    start=0;
    for(n=0; n<numField; n++)
      for(i=1; i<share; i++)
        for(j=0; j<ny; j++) 
        {
          for(k=0; k<nz; k++) {
	  	      data[start+k]=D->shareF[n][istart-i][j][k]; 
            D->shareF[n][istart-i][j][k]=0.0; 
          }
	  	    start+=nz; 
	  	  }    
    MPI_Send(data,num,MPI_DOUBLE,D->prevXrank,myrank, MPI_COMM_WORLD);
  }
  MPI_Barrier(MPI_COMM_WORLD);
  free(data);
}
                                                
                                        
void MPI_TransferJ_Period_X(Domain *D,int numField)
{
  int n,i,j,k,start,ny,nz,share;
  int istart,iend,jstart,jend,kstart,kend;
  int myrank, nTasks, rank,num;
  double *data;

  MPI_Status status;

  MPI_Comm_size(MPI_COMM_WORLD, &nTasks);
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

  istart=D->istart;    iend=D->iend;
  jstart=D->jstart;    jend=D->jend;
  kstart=D->kstart;    kend=D->kend;
	ny=1;	nz=1;
	if(D->dimension>1) {	ny=jend+3; } else ;
	if(D->dimension>2) {	nz=kend+3; } else ;

  rank=myrank/(D->M*D->N);
  share=D->shareX;

	// Right to Left 
  num=numField*ny*nz*share;
  data = (double *)malloc(num*sizeof(double ));

  start=0;
  for(n=0; n<numField; n++)
    for(i=0; i<share; i++)
      for(j=0; j<ny; j++)  
      {
        for(k=0; k<nz; k++)  {
          data[start+k]=D->shareF[n][iend+i][j][k]; 
          D->shareF[n][iend+i][j][k]=0.0; 
        }
        start+=nz; 
      }

	if(D->L==1) {
    start=0;
    for(n=0; n<numField; n++)
      for(i=0; i<share; i++)
        for(j=0; j<ny; j++) {
          for(k=0; k<nz; k++) 
            D->shareF[n][istart+i][j][k]+=data[k+start]; 
          start+=nz; 
        }
			
	} else {
		if(rank==0)   {
      MPI_Recv(data,num,MPI_DOUBLE,myrank+(D->L-1)*D->M*D->N,myrank+(D->L-1)*D->M*D->N,MPI_COMM_WORLD,&status);
      start=0;
      for(n=0; n<numField; n++)
        for(i=0; i<share; i++)
          for(j=0; j<ny; j++) {
            for(k=0; k<nz; k++) 
              D->shareF[n][istart+i][j][k]+=data[k+start]; 
            start+=nz; 
          }
    }  else if(rank==D->L-1)
      MPI_Send(data,num,MPI_DOUBLE,myrank-(D->L-1)*D->M*D->N,myrank, MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);
	}
	free(data);

	// Left to Right 
  num=numField*ny*nz*2;
  data = (double *)malloc(num*sizeof(double ));

  start=0;
  for(n=0; n<numField; n++)
    for(i=1; i<share; i++)
      for(j=0; j<ny; j++) 
      {
        for(k=0; k<nz; k++) {
          data[start+k]=D->shareF[n][istart-i][j][k]; 
          D->shareF[n][istart-i][j][k]=0.0;
        } 
        start+=nz; 
      }

	if(D->L==1) {
    start=0;
    for(n=0; n<numField; n++)
      for(i=1; i<share; i++)
        for(j=0; j<ny; j++) {
          for(k=0; k<nz; k++) 
            D->shareF[n][iend-i][j][k]+=data[k+start];
          start+=nz; 
        }
	} else {
		if(rank==D->L-1) {
      MPI_Recv(data,num,MPI_DOUBLE,myrank-(D->L-1)*D->M*D->N,myrank-(D->L-1)*D->M*D->N,MPI_COMM_WORLD,&status);
      start=0;
      for(n=0; n<numField; n++)
        for(i=1; i<share; i++)
          for(j=0; j<ny; j++) {
            for(k=0; k<nz; k++) 
              D->shareF[n][iend-i][j][k]+=data[k+start]; 
            start+=nz; 
          }
    } else if(rank==0)
      MPI_Send(data,num,MPI_DOUBLE,myrank+(D->L-1)*D->M*D->N,myrank, MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);
	}
   
  free(data);
}
